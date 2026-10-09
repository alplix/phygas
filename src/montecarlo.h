// Monte Carlo resampling for orbital/astrometric uncertainty propagation.
// Scheme decided 2026-10-06 by reading Bailer-Jones et al. 2019
// (arXiv:1911.02473, the actual precedent this project scales up), not
// guessed -- see NOTES_method_selection.md's "Monte Carlo sampling scheme"
// section for the full reasoning and the two-stage design this feeds into.
//
// Deliberately a fully deterministic, seedable PRNG (not libc rand(), which
// has no specified quality/period and can't be guaranteed reproducible
// across platforms) -- every BOINC work unit needs its Monte Carlo draws to
// be exactly reproducible from the same seed for result validation (quorum
// agreement between independently-run replicas depends on this).
#ifndef NOSTOS_MONTECARLO_H
#define NOSTOS_MONTECARLO_H

#include <stdint.h>

// splitmix64 (Vigna) -- the PRNG itself. Tiny, fast, passes standard
// statistical test suites (BigCrush) despite its simplicity; used here
// directly rather than as just a seed-expander for something fancier,
// since this project's draw volume per work unit (hundreds of Gaussian
// samples, not billions) doesn't need a larger-period generator.
typedef struct {
    uint64_t state;
} McRng;

static inline void mc_rng_seed(McRng *rng, uint64_t seed) {
    rng->state = seed;
}

static inline uint64_t mc_rng_next_u64(McRng *rng) {
    uint64_t z = (rng->state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

// Uniform double in [0, 1). 53 bits of the 64-bit output -> exact double
// precision, no bias from a naive (double)rand()/RAND_MAX-style division.
static inline double mc_rng_uniform(McRng *rng) {
    return (double)(mc_rng_next_u64(rng) >> 11) * (1.0 / 9007199254740992.0); // 2^53
}

// Standard normal (mean 0, variance 1) via Box-Muller. Returns one value
// per call (the "spare" second value Box-Muller naturally produces is
// discarded rather than cached -- simpler and one call = one independent
// draw, deliberately not optimizing away half the trig calls here since
// draw volume is small; see this header's top comment).
double mc_gaussian(McRng *rng);

// Maximum dimension this module's multivariate sampler supports -- fixed-
// size, no dynamic allocation, matching this project's style elsewhere
// (Vec3 etc.). 6 covers the largest real use case (a star's 5 Gaia
// astrometric parameters -- parallax, pmra, pmdec, ra, dec -- plus radial
// velocity), with headroom.
#define MC_MAX_DIM 6

// Cholesky decomposition of a symmetric positive-definite covariance matrix
// cov[n][n] (only the lower triangle is read) into L such that L*L^T = cov.
// Returns 0 on success, nonzero if cov isn't positive-definite (a zero or
// negative value appears on the diagonal during decomposition -- callers
// should treat this as "bad input", e.g. a mis-specified covariance matrix,
// not silently proceed).
int mc_cholesky(const double cov[MC_MAX_DIM][MC_MAX_DIM], int n, double L[MC_MAX_DIM][MC_MAX_DIM]);

// Draws one correlated multivariate-Gaussian sample: out = mean + L*z,
// where z is a vector of n independent standard normals and L is the
// Cholesky factor from mc_cholesky() (so out has the originally-specified
// mean and covariance). This is the standard technique for sampling any
// multivariate Gaussian from independent univariate ones -- used here for
// both an ISO's 3x3 velocity-uncertainty covariance (Bailer-Jones' "clones")
// and a star's correlated Gaia astrometric errors, with the SAME function:
// the two use cases differ only in which covariance matrix is passed in,
// not in any sampling logic.
void mc_sample_multivariate_gaussian(
    McRng *rng, const double mean[MC_MAX_DIM], const double L[MC_MAX_DIM][MC_MAX_DIM],
    int n, double out[MC_MAX_DIM]
);

#endif
