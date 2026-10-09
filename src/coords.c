// See coords.h for the public interface and its caveat about independent validation still being
// needed. This file documents every constant's source so a future cross-check knows exactly
// what to verify.
//
// Galactocentric frame convention used throughout this file (our own choice, not universal --
// matches the shape of the common astropy/galpy convention but documented here independently so
// internal consistency doesn't depend on getting someone else's convention exactly right from
// memory):
//   - Right-handed Cartesian, Galactic center at the origin.
//   - Sun at (-r0, 0, z0).
//   - +X points from the Sun towards the Galactic center.
//   - +Y points in the direction of Galactic rotation at the Sun's location.
//   - +Z points towards the North Galactic Pole.
#include "coords.h"
#include "units.h"
#include <math.h>

#define DEG2RAD (3.14159265358979323846 / 180.0)
#define RAD2DEG (180.0 / 3.14159265358979323846)

// ICRS -> Galactic Cartesian rotation matrix. This exact 3x3 orthogonal matrix (derived from the
// Hipparcos-era definition of the Galactic pole/node in ICRS: alpha_NGP=192.85948deg,
// delta_NGP=27.12825deg, l_NCP=122.93192deg) is one of the most widely reproduced constants in
// astronomical software (ESA 1997 Hipparcos & Tycho Catalogues Vol 1 Sec 1.5.3; reproduced
// verbatim by astropy and many other packages).
//
// Cross-checked 2026-10-06: tests/test_horizons_crosscheck.c applies this matrix's transpose to
// orbital_elements.c's output and recovers a real JPL Horizons ICRS state vector for 1I/'Oumuamua
// to ~1e-6 relative error -- confirming both the digits AND the sign/orientation convention are
// right, not just "plausibly a rotation matrix".
static const double A_G[3][3] = {
    {-0.0548755604162154, -0.8734370902348850, -0.4838350155487132},
    { 0.4941094278755837, -0.4448296299600112,  0.7469822444972189},
    {-0.8676661490190047, -0.1980763734312015,  0.4559837761750669},
};

static Vec3 mat_vec(const double m[3][3], Vec3 v) {
    return (Vec3){
        m[0][0]*v.x + m[0][1]*v.y + m[0][2]*v.z,
        m[1][0]*v.x + m[1][1]*v.y + m[1][2]*v.z,
        m[2][0]*v.x + m[2][1]*v.y + m[2][2]*v.z,
    };
}

static Vec3 mat_transpose_vec(const double m[3][3], Vec3 v) {
    return (Vec3){
        m[0][0]*v.x + m[1][0]*v.y + m[2][0]*v.z,
        m[0][1]*v.x + m[1][1]*v.y + m[2][1]*v.z,
        m[0][2]*v.x + m[1][2]*v.y + m[2][2]*v.z,
    };
}

// Sun's assumed position/velocity. r0: GRAVITY Collaboration 2018 (Sgr A* orbit). v_circ,
// u/v/w_pec: Schonrich, Binney & Dehnen 2010 (their best-fit LSR circular speed and solar
// peculiar motion, U positive towards the Galactic center). z0: a commonly cited ~tens-of-pc
// figure, less critical than r0/v_circ for this project's purposes. Like the Milky Way
// potential in potential.c, these are standard literature values, not yet cross-checked here.
const SolarParams solar_params_default = {
    .r0_kpc = 8.122,
    .z0_kpc = 0.0208,
    .v_circ_kms = 229.0,
    .u_pec_kms = 11.1,
    .v_pec_kms = 12.24,
    .w_pec_kms = 7.25,
};

// AU/yr -> km/s, derived from first principles (not a recalled external constant): 1 AU/yr =
// 1.495978707e8 km / 3.15576e7 s (Julian year). By the literal definition of the parsec, a
// proper motion of 1 arcsec/yr at a distance of 1 pc corresponds to a transverse velocity of
// exactly 1 AU/yr -- so this same constant converts (arcsec/yr * pc) directly to km/s.
#define AUYR_TO_KMS (1.495978707e8 / 3.15576e7)

static Vec3 sun_velocity_galcen(const SolarParams *sun) {
    return (Vec3){ sun->u_pec_kms, sun->v_circ_kms + sun->v_pec_kms, sun->w_pec_kms };
}

void astrometry_to_galactocentric(Astrometry a, const SolarParams *sun, PhaseSpace *out) {
    double ra = a.ra_deg * DEG2RAD;
    double dec = a.dec_deg * DEG2RAD;
    double d_kpc = 1.0 / a.parallax_mas; // see header: mas parallax -> kpc distance is a clean 1/x
    double d_pc = d_kpc * 1000.0;

    Vec3 r_hat = { cos(dec)*cos(ra), cos(dec)*sin(ra), sin(dec) };
    Vec3 ra_hat = { -sin(ra), cos(ra), 0.0 };
    Vec3 dec_hat = { -sin(dec)*cos(ra), -sin(dec)*sin(ra), cos(dec) };

    // Position: scale the ICRS-aligned unit vector by distance (still heliocentric, not yet
    // rotated to Galactic or shifted to Galactocentric -- that part is shared with
    // orbital_elements.c via helio_equatorial_to_galactocentric_pos()).
    Vec3 pos_icrs_kpc = v3_scale(r_hat, d_kpc);
    out->pos = helio_equatorial_to_galactocentric_pos(pos_icrs_kpc, sun);

    // Velocity: proper motion (mas/yr) + radial velocity (km/s) -> Cartesian kpc/Myr in the
    // ICRS-aligned basis, then the same shared rotate+shift helper as position.
    double v_ra_kms = AUYR_TO_KMS * 1e-3 * a.pmra_masyr * d_pc;
    double v_dec_kms = AUYR_TO_KMS * 1e-3 * a.pmdec_masyr * d_pc;
    Vec3 v_icrs_kms = v3_add(v3_add(v3_scale(r_hat, a.rv_kms), v3_scale(ra_hat, v_ra_kms)), v3_scale(dec_hat, v_dec_kms));
    Vec3 v_icrs_kpcmyr = {
        nostos_kms_to_kpcmyr(v_icrs_kms.x),
        nostos_kms_to_kpcmyr(v_icrs_kms.y),
        nostos_kms_to_kpcmyr(v_icrs_kms.z),
    };
    out->vel = helio_equatorial_to_galactocentric_vel(v_icrs_kpcmyr, sun);
}

Vec3 helio_equatorial_to_galactocentric_pos(Vec3 pos_icrs_kpc, const SolarParams *sun) {
    Vec3 helio_gal = mat_vec(A_G, pos_icrs_kpc);
    return (Vec3){
        helio_gal.x - sun->r0_kpc,
        helio_gal.y,
        helio_gal.z + sun->z0_kpc,
    };
}

Vec3 helio_equatorial_to_galactocentric_vel(Vec3 vel_icrs_kpcmyr, const SolarParams *sun) {
    Vec3 v_gal_helio = mat_vec(A_G, vel_icrs_kpcmyr);
    Vec3 v_sun_kms = sun_velocity_galcen(sun);
    Vec3 v_sun_kpcmyr = {
        nostos_kms_to_kpcmyr(v_sun_kms.x),
        nostos_kms_to_kpcmyr(v_sun_kms.y),
        nostos_kms_to_kpcmyr(v_sun_kms.z),
    };
    return v3_add(v_gal_helio, v_sun_kpcmyr);
}

Astrometry galactocentric_to_astrometry(PhaseSpace ps, const SolarParams *sun) {
    Vec3 helio_gal = {
        ps.pos.x + sun->r0_kpc,
        ps.pos.y,
        ps.pos.z - sun->z0_kpc,
    };
    double d_kpc = v3_norm(helio_gal);
    double d_pc = d_kpc * 1000.0;
    Vec3 p_gal_unit = v3_scale(helio_gal, 1.0 / d_kpc);

    Vec3 r_hat = mat_transpose_vec(A_G, p_gal_unit); // ICRS unit vector
    double ra = atan2(r_hat.y, r_hat.x);
    double dec = asin(r_hat.z);
    if (ra < 0) ra += 2.0 * 3.14159265358979323846;

    Vec3 ra_hat = { -sin(ra), cos(ra), 0.0 };
    Vec3 dec_hat = { -sin(dec)*cos(ra), -sin(dec)*sin(ra), cos(dec) };

    Vec3 v_sun = sun_velocity_galcen(sun);
    Vec3 v_galcen_kms = {
        nostos_kpcmyr_to_kms(ps.vel.x),
        nostos_kpcmyr_to_kms(ps.vel.y),
        nostos_kpcmyr_to_kms(ps.vel.z),
    };
    Vec3 v_gal_helio = v3_sub(v_galcen_kms, v_sun);
    Vec3 v_icrs = mat_transpose_vec(A_G, v_gal_helio);

    Astrometry a;
    a.ra_deg = ra * RAD2DEG;
    a.dec_deg = dec * RAD2DEG;
    a.parallax_mas = 1.0 / d_kpc;
    a.rv_kms = v3_dot(v_icrs, r_hat);
    a.pmra_masyr = v3_dot(v_icrs, ra_hat) / (AUYR_TO_KMS * 1e-3 * d_pc);
    a.pmdec_masyr = v3_dot(v_icrs, dec_hat) / (AUYR_TO_KMS * 1e-3 * d_pc);
    return a;
}
