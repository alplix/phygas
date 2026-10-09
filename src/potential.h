// Generic Galactic-potential interface: an acceleration law (per unit mass, so
// it applies equally to a Gaia star or an interstellar object -- orbits in a
// smooth background potential don't depend on the test particle's own mass)
// plus a scalar potential energy, used for energy-conservation checks.
//
// All physical quantities are in the internal unit system from units.h:
// kpc, Myr, kpc/Myr, Msun.
#ifndef NOSTOS_POTENTIAL_H
#define NOSTOS_POTENTIAL_H

#include "vec3.h"

typedef struct Potential Potential;

struct Potential {
    void (*acceleration)(const Potential *pot, Vec3 pos, Vec3 *accel_out);
    double (*energy)(const Potential *pot, Vec3 pos);
    const void *params;
};

static inline Vec3 potential_acceleration(const Potential *pot, Vec3 pos) {
    Vec3 a;
    pot->acceleration(pot, pos, &a);
    return a;
}

static inline double potential_energy(const Potential *pot, Vec3 pos) {
    return pot->energy(pot, pos);
}

// -- Point-mass (Kepler) potential --------------------------------------
// Phi(r) = -GM/r. Exists purely to validate the integrator itself (closed
// orbits, exact conservation, exact analytic position-vs-time via Kepler's
// equation) completely independent of any Milky-Way parameter uncertainty.
typedef struct {
    double GM; // kpc^3 / Myr^2
} KeplerParams;

Potential kepler_potential(const KeplerParams *params);

// -- Milky Way potential: Allen & Santillan (1991) --------------------
// Replaced the earlier generic Hernquist+Miyamoto-Nagai+NFW placeholder
// (2026-10-06) with Allen & Santillan's actual published model (A&SS 86,
// 1991, "An improved model of the galactic mass distribution for orbit
// computations"), a widely-used closed-form 3-component model -- the same
// one used for prior stellar-flyby/orbit-backtracking work this project's
// own NOTES_method_selection.md cites (e.g. Garcia-Sanchez et al. 2001).
//
// Three components, confirmed against multiple independent secondary
// sources (mass unit M_gal = 2.325e7 Msun; parameter table values 606.0 /
// 3690.0 / 4615.0 M_gal and 0.3873 / 5.3178 / 0.2500 / 12.0 kpc; halo cutoff
// Lambda=100 kpc, gamma=2.02 -- see potential.c's comment for the derivation
// of the halo's closed-form potential from these):
//   - Bulge: PLUMMER, not Hernquist -- a real shape difference from the
//     earlier placeholder, not just new numbers. Phi = -GM/sqrt(r^2+b^2).
//   - Disk: Miyamoto-Nagai (same functional form the placeholder already
//     used, just new mass/scale values).
//   - Halo: a generalized truncated enclosed-mass law, M(r) =
//     M_h*(r/a_h)^gamma / (1+(r/a_h)^(gamma-1)) for r<Lambda, M(Lambda)
//     (constant) beyond the cutoff -- reduces to a plain point mass outside
//     Lambda by Newton's shell theorem, which is how the acceleration is
//     computed either way (same -GM(r)/r^2*r_hat pattern the old NFW used).
//
// tests/test_milkyway_allen_santillan.c checks: (a) the halo's energy and
// acceleration are self-consistent (numerical dPhi/dr matches the analytic
// acceleration) across the r=Lambda truncation boundary, where a sign or
// continuity mistake would most likely show up; (b) the resulting rotation
// curve falls in the ~220-235 km/s range the model is known for near the
// solar radius. The disk's Miyamoto-Nagai formula itself is still the one
// already cross-checked against galpy in tests/test_milkyway_galpy.c (that
// test's bulge/halo checks no longer apply now that those shapes changed --
// see its own updated comment).
typedef struct {
    double bulge_GM;    // kpc^3/Myr^2
    double bulge_b;     // kpc, Plummer scale radius
    double disk_GM;     // kpc^3/Myr^2
    double disk_a;      // kpc, Miyamoto-Nagai radial scale length
    double disk_b;      // kpc, Miyamoto-Nagai vertical scale height
    double halo_GM;     // kpc^3/Myr^2, G*M_h (the halo's characteristic mass scale)
    double halo_a;      // kpc, halo scale radius
    double halo_cutoff; // kpc, truncation radius (Lambda)
    double halo_gamma;  // dimensionless enclosed-mass exponent (A&S: 2.02)
} MilkyWayParams;

extern const MilkyWayParams milkyway_potential_default;

Potential milkyway_potential(const MilkyWayParams *params);

#endif
