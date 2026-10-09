#include "potential.h"
#include "units.h"
#include <math.h>

// ---- Kepler -------------------------------------------------------------

static void kepler_accel(const Potential *pot, Vec3 pos, Vec3 *out) {
    const KeplerParams *p = (const KeplerParams *)pot->params;
    double r = v3_norm(pos);
    double f = -p->GM / (r * r * r); // a = -GM/r^3 * r_vec
    *out = v3_scale(pos, f);
}

static double kepler_energy(const Potential *pot, Vec3 pos) {
    const KeplerParams *p = (const KeplerParams *)pot->params;
    double r = v3_norm(pos);
    return -p->GM / r;
}

Potential kepler_potential(const KeplerParams *params) {
    Potential pot;
    pot.acceleration = kepler_accel;
    pot.energy = kepler_energy;
    pot.params = params;
    return pot;
}

// ---- Milky Way: Allen & Santillan (1991) --------------------------------
// See potential.h's big comment for the model's provenance and the
// parameter table's sources. *_GM fields are G_internal * mass_in_Msun
// (mass unit M_gal = 2.325e7 Msun, per A&S's own convention) -- computed
// here rather than left as raw masses so callers never need to know the
// unit-conversion constant.
#define NOSTOS_MGAL_TO_MSUN 2.325e7
const MilkyWayParams milkyway_potential_default = {
    .bulge_GM = NOSTOS_G * (606.0 * NOSTOS_MGAL_TO_MSUN),
    .bulge_b = 0.3873,
    .disk_GM = NOSTOS_G * (3690.0 * NOSTOS_MGAL_TO_MSUN),
    .disk_a = 5.3178,
    .disk_b = 0.2500,
    .halo_GM = NOSTOS_G * (4615.0 * NOSTOS_MGAL_TO_MSUN),
    .halo_a = 12.0,
    .halo_cutoff = 100.0,
    .halo_gamma = 2.02,
};

// Halo enclosed mass (as G*M, i.e. already in this project's *_GM
// convention), A&S's generalized truncated profile: for r<Lambda it's a
// generalized-NFW-like power law; beyond the cutoff the enclosed mass is
// frozen at its value at Lambda (Newton's shell theorem then makes the
// exterior potential/acceleration identical to a point mass of that size --
// this single function covers both branches, used by both milkyway_accel
// and milkyway_energy so they can never disagree about where the cutoff is).
static double halo_enclosed_GM(const MilkyWayParams *p, double r) {
    double r_eff = (r < p->halo_cutoff) ? r : p->halo_cutoff;
    double x = r_eff / p->halo_a;
    return p->halo_GM * pow(x, p->halo_gamma) / (1.0 + pow(x, p->halo_gamma - 1.0));
}

static void milkyway_accel(const Potential *pot, Vec3 pos, Vec3 *out) {
    const MilkyWayParams *p = (const MilkyWayParams *)pot->params;

    double R2 = pos.x * pos.x + pos.y * pos.y;
    double R = sqrt(R2);
    double r = v3_norm(pos);

    Vec3 a_total = {0, 0, 0};

    // Plummer bulge: Phi = -GM/sqrt(r^2+b^2); accel = -GM/(r^2+b^2)^1.5 * pos
    {
        double denom = pow(r * r + p->bulge_b * p->bulge_b, 1.5);
        double f = -p->bulge_GM / (denom > 0 ? denom : 1e-12);
        a_total = v3_add(a_total, v3_scale(pos, f));
    }

    // Miyamoto-Nagai disk (same formula as before, new parameter values).
    {
        double zeta = sqrt(pos.z * pos.z + p->disk_b * p->disk_b);
        double azeta = p->disk_a + zeta;
        double D2 = R2 + azeta * azeta;
        double D = sqrt(D2);
        double D3 = D2 * D;
        double aR = -p->disk_GM * R / D3; // acceleration in the R direction
        double az = -p->disk_GM * pos.z * azeta / (zeta > 0 ? zeta : 1e-12) / D3;
        if (R > 1e-12) {
            a_total.x += aR * pos.x / R;
            a_total.y += aR * pos.y / R;
        }
        a_total.z += az;
    }

    // Halo: accel = -G*M_enc(r)/r^2 * r_hat, same shell-theorem pattern the
    // earlier NFW placeholder used, just with A&S's enclosed-mass law.
    {
        double r_safe = (r > 1e-12 ? r : 1e-12);
        double f = -halo_enclosed_GM(p, r_safe) / (r_safe * r_safe * r_safe);
        a_total = v3_add(a_total, v3_scale(pos, f));
    }

    *out = a_total;
}

static double milkyway_energy(const Potential *pot, Vec3 pos) {
    const MilkyWayParams *p = (const MilkyWayParams *)pot->params;

    double R2 = pos.x * pos.x + pos.y * pos.y;
    double r = v3_norm(pos);

    double phi_bulge = -p->bulge_GM / sqrt(r * r + p->bulge_b * p->bulge_b);

    double zeta = sqrt(pos.z * pos.z + p->disk_b * p->disk_b);
    double azeta = p->disk_a + zeta;
    double D = sqrt(R2 + azeta * azeta);
    double phi_disk = -p->disk_GM / D;

    // Halo: Phi(r) = -GM_enc(Lambda)/Lambda for r>=Lambda (point-mass
    // exterior, by the shell theorem); for r<Lambda, integrate the
    // shell-theorem force -GM_enc(s)/s^2 from r out to Lambda and add that
    // exterior value (derivation: substitute x=s/a_h in the generalized
    // enclosed-mass law, giving a clean log antiderivative -- see
    // NOTES_method_selection.md for the full derivation written out).
    double r_safe = (r > 1e-12 ? r : 1e-12);
    double GM_Lambda = halo_enclosed_GM(p, p->halo_cutoff); // = halo_GM at r=Lambda exactly
    double phi_halo;
    if (r_safe >= p->halo_cutoff) {
        phi_halo = -GM_Lambda / r_safe;
    } else {
        double gm1 = p->halo_gamma - 1.0;
        double xL = p->halo_cutoff / p->halo_a;
        double xr = r_safe / p->halo_a;
        double integral = (1.0 / gm1) * (log(1.0 + pow(xL, gm1)) - log(1.0 + pow(xr, gm1)));
        phi_halo = -(p->halo_GM / p->halo_a) * integral - GM_Lambda / p->halo_cutoff;
    }

    return phi_bulge + phi_disk + phi_halo;
}

Potential milkyway_potential(const MilkyWayParams *params) {
    Potential pot;
    pot.acceleration = milkyway_accel;
    pot.energy = milkyway_energy;
    pot.params = params;
    return pot;
}
