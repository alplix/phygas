// Internal unit system for the whole integrator: length in kpc, time in Myr,
// velocity in kpc/Myr, mass in solar masses. Fully self-consistent (no hidden
// factors inside force/integrator code) -- convert at the edges only, with the
// helpers below, when data comes in as km/s (the usual astronomical unit).
#ifndef NOSTOS_UNITS_H
#define NOSTOS_UNITS_H

// 1 kpc / 1 Myr, expressed in km/s. kpc->km and Myr->s are both exact SI
// definitions (Julian year = 365.25 days), so this is exact to the precision
// shown, not a fitted/rounded value.
//   1 kpc = 3.0856775814913673e16 km
//   1 Myr = 3.15576e13 s   (365.25 * 86400 * 1e6)
//   (km/s) per (kpc/Myr) = 3.0856775814913673e16 / 3.15576e13
#define NOSTOS_KMS_PER_KPCMYR 977.7922216807892

static inline double nostos_kms_to_kpcmyr(double v_kms) {
    return v_kms / NOSTOS_KMS_PER_KPCMYR;
}
static inline double nostos_kpcmyr_to_kms(double v_kpcmyr) {
    return v_kpcmyr * NOSTOS_KMS_PER_KPCMYR;
}

// Newton's constant, G = 4.30091e-6 kpc*(km/s)^2 / Msun (standard value quoted
// e.g. by galpy/astropy; derivable from G = 6.6743e-11 m^3 kg^-1 s^-2 plus
// Msun = 1.98892e30 kg and 1 kpc = 3.0856775814913673e19 m).
// Converted here into the internal unit system: kpc^3 Msun^-1 Myr^-2.
#define NOSTOS_G_KPC_KMS2_PER_MSUN 4.30091e-6
#define NOSTOS_G (NOSTOS_G_KPC_KMS2_PER_MSUN / (NOSTOS_KMS_PER_KPCMYR * NOSTOS_KMS_PER_KPCMYR))

// 1 AU in kpc. Exact by definition: 1 pc = 206264.80624709636 AU (a parsec is defined as the
// distance at which 1 AU subtends exactly 1 arcsec, and 1 radian = 206264.80624709636 arcsec).
#define NOSTOS_AU_PER_KPC (206264.80624709636 * 1000.0)
#define NOSTOS_KPC_PER_AU (1.0 / NOSTOS_AU_PER_KPC)

// 1 day in Myr, consistent with this project's Julian-year convention (365.25 days/yr) --
// orbital-element epochs (time of perihelion passage etc.) are conventionally given in days.
#define NOSTOS_MYR_PER_DAY (1.0 / (365.25 * 1.0e6))

// GM of the Sun in this unit system: the Sun's mass IS 1 Msun by definition, so this is just
// NOSTOS_G itself -- no separate GM_sun constant to keep in sync.
#define NOSTOS_GM_SUN NOSTOS_G

// Total GM of the solar system (Sun + planets + ...) relative to the Sun's: the "Keplerian GM" JPL
// Horizons uses for orbital elements about the Solar System Barycentre (DE441), 2.9630927493968080e-4
// au^3/d^2, divided by the Gaussian k^2 = 2.9591220828559115e-4 au^3/d^2. Needed with barycentric
// elements (see orbital_elements_inbound_asymptote_gm()).
#define NOSTOS_GM_SSB (NOSTOS_GM_SUN * (2.9630927493968080e-4 / 2.9591220828559115e-4))

#endif
