#include "montecarlo.h"
#include <math.h>

double mc_gaussian(McRng *rng) {
    // Box-Muller (basic form, not the polar/Marsaglia variant -- no
    // rejection loop needed, and this project's draw volume is small
    // enough that the extra sin/cos call isn't worth the added complexity
    // of caching the second value).
    double u1 = mc_rng_uniform(rng);
    // u1 == 0 would make log(u1) = -inf; vanishingly rare (probability
    // ~2^-53) but cheap to guard against outright rather than rely on that
    // rarity.
    if (u1 <= 0.0) u1 = 1.0 / 18446744073709551616.0; // smallest representable nonzero draw
    double u2 = mc_rng_uniform(rng);
    double r = sqrt(-2.0 * log(u1));
    return r * cos(2.0 * 3.14159265358979323846 * u2);
}

int mc_cholesky(const double cov[MC_MAX_DIM][MC_MAX_DIM], int n, double L[MC_MAX_DIM][MC_MAX_DIM]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = cov[i][j];
            for (int k = 0; k < j; k++) {
                sum -= L[i][k] * L[j][k];
            }
            if (i == j) {
                if (sum <= 0.0) return 1; // not positive-definite
                L[i][j] = sqrt(sum);
            } else {
                L[i][j] = sum / L[j][j];
            }
        }
        for (int j = i + 1; j < n; j++) {
            L[i][j] = 0.0; // upper triangle unused, but zero it so callers
                            // can't accidentally read garbage if n < MC_MAX_DIM
        }
    }
    return 0;
}

void mc_sample_multivariate_gaussian(
    McRng *rng, const double mean[MC_MAX_DIM], const double L[MC_MAX_DIM][MC_MAX_DIM],
    int n, double out[MC_MAX_DIM]
) {
    double z[MC_MAX_DIM];
    for (int i = 0; i < n; i++) {
        z[i] = mc_gaussian(rng);
    }
    for (int i = 0; i < n; i++) {
        double sum = mean[i];
        for (int j = 0; j <= i; j++) {
            sum += L[i][j] * z[j];
        }
        out[i] = sum;
    }
}
