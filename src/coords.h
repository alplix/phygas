// Observed astrometry (ICRS RA/Dec/parallax/proper motion/radial velocity -- what Gaia and the
// ISO orbital-element solutions actually give us) <-> Galactocentric Cartesian phase space (what
// the integrator in integrator.h/potential.h needs).
//
// The round-trip self-consistency test in tests/test_coords.c (forward transform then exact
// inverse, recovering the original input to float precision) catches a sign error or a
// non-orthogonal matrix, but can't catch an internally-consistent-but-wrong convention. Two
// independent checks now exist: tests/test_horizons_crosscheck.c runs this module's shared
// helio_equatorial_to_galactocentric_{pos,vel}() helpers (via orbital_elements.c) against a real
// JPL Horizons state vector for 1I/'Oumuamua (~1e-6 relative error), and
// tests/test_coords_astropy_crosscheck.c runs astrometry_to_galactocentric() ITSELF (the
// star-side function, not shared with orbital_elements.c) against astropy -- a real, independent,
// widely-used implementation of the same ICRS<->Galactic transform -- for a real star (HIP 3757),
// matching to ~1e-7 relative error. See coords.c's own comment on A_G. The Milky Way potential in
// potential.c is a separate subsystem with its own cross-check (test_milkyway_galpy.c/
// test_milkyway_allen_santillan.c).
#ifndef NOSTOS_COORDS_H
#define NOSTOS_COORDS_H

#include "vec3.h"

typedef struct {
    double ra_deg;       // ICRS right ascension, degrees
    double dec_deg;      // ICRS declination, degrees
    double parallax_mas; // parallax, milliarcseconds (must be > 0)
    double pmra_masyr;   // proper motion in RA, mas/yr (already *cos(dec), Gaia's own convention)
    double pmdec_masyr;  // proper motion in Dec, mas/yr
    double rv_kms;       // radial velocity, km/s (positive = receding)
} Astrometry;

typedef struct {
    Vec3 pos; // kpc
    Vec3 vel; // kpc/Myr
} PhaseSpace;

// Sun's assumed position and velocity in the Galactocentric frame this module defines (see
// coords.c's top comment for the exact frame convention: right-handed, Galactic center at
// origin, Sun at (-r0, 0, z0), rotation in the +y direction at the Sun's location).
// Defaults are commonly-cited modern values (GRAVITY Collaboration 2018 for r0; Schonrich et al.
// 2010 for the Sun's peculiar motion) -- like the Milky Way potential, NOT yet cross-checked
// against a reference implementation, just standard literature values.
typedef struct {
    double r0_kpc;     // Sun's Galactocentric radius
    double z0_kpc;      // Sun's height above the Galactic plane
    double v_circ_kms;  // circular rotation speed at r0
    double u_pec_kms;   // Sun's peculiar motion relative to the LSR, towards the Galactic center
    double v_pec_kms;   // ... in the direction of Galactic rotation
    double w_pec_kms;   // ... towards the North Galactic Pole
} SolarParams;

extern const SolarParams solar_params_default;

// Forward: observed astrometry -> Galactocentric phase space.
void astrometry_to_galactocentric(Astrometry a, const SolarParams *sun, PhaseSpace *out);

// Inverse: Galactocentric phase space -> observed astrometry (as seen from the Sun's assumed
// position/velocity). Exists mainly so tests/test_coords.c can round-trip a made-up Astrometry
// through both functions and check it comes back unchanged -- the main self-consistency check
// this module has.
Astrometry galactocentric_to_astrometry(PhaseSpace ps, const SolarParams *sun);

// Shared building blocks for any OTHER heliocentric-equatorial-Cartesian source (e.g.
// orbital_elements.c, which gets a solar-system body's position from its orbital elements
// rather than from astrometry) to reach the same Galactocentric frame without duplicating the
// rotate-then-shift logic. pos/vel are heliocentric, ICRS-axis-aligned (NOT yet rotated to
// Galactic), in kpc / kpc-per-Myr.
Vec3 helio_equatorial_to_galactocentric_pos(Vec3 pos_icrs_kpc, const SolarParams *sun);
Vec3 helio_equatorial_to_galactocentric_vel(Vec3 vel_icrs_kpcmyr, const SolarParams *sun);

#endif
