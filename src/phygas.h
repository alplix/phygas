// Phygas: where did runaway stars and neutron stars come from?
//
// Every target (a Gaia star or a pulsar) is traced back through the Milky Way together with every
// young star cluster it could have reached, both drawn many times from their measurement
// uncertainties (Monte Carlo on both sides), and the closest approach of each pair is recorded.
// A target born in a cluster passes through it at a time no larger than the cluster's age.
//
// Units inside: kpc, Myr, kpc/Myr (as in the Nostos core this reuses: coords/potential/montecarlo).
// Reported: pc, Myr (positive = in the past), km/s.
#ifndef PHYGAS_H
#define PHYGAS_H

#include <stdint.h>
#include "coords.h"
#include "montecarlo.h"
#include "potential.h"
#include "vec3.h"

typedef struct {
    char name[40];
    char kind;              // o = open cluster, m = moving group
    double ra, dec;         // deg
    double dist_pc, e_dist_pc;
    double pmra, e_pmra, pmdec, e_pmdec;   // mas/yr (pmra includes cos dec)
    double rv, e_rv;        // km/s
    int has_rv;             // 0: drawn from the Galactic-rotation expectation +- PH_CLUSTER_RV_PRIOR
    double r50_pc, rtot_pc; // radius holding half the members, total radius
    double age16, age50, age84;  // Myr
    double tmax_myr;        // how far back to trace for this cluster
} PhCluster;

typedef struct {
    char id[40];
    char kind;              // P = pulsar / neutron star, S = star
    double ra, dec;
    double plx, e_plx;      // mas; used when plx > 0
    double dist_pc, e_lndist;  // used when plx == 0: log-normal distance
    double pmra, e_pmra, pmdec, e_pmdec;
    double rv, e_rv;
    int has_rv;             // pulsars: never; drawn around the rotation expectation +- PH_PULSAR_RV_PRIOR
} PhTarget;

// Nominal path and its uncertainty, for the prefilter.
typedef struct {
    int n;                  // points, PH_PRE_STEP_MYR apart, from now into the past
    Vec3 *pos;              // kpc, Galactocentric
    double sig_pos;         // kpc, spread of the present position across clones
    double sig_vel;         // kpc/Myr, spread of the velocity across clones
} PhTrack;

typedef struct {
    int n;                  // clones run
    double d05, d50, d95;   // closest approach, pc
    double t05, t50, t95;   // its time, Myr ago
    double v05, v50, v95;   // relative speed at closest approach, km/s (the ejection velocity)
    double f_r50, f_rtot;   // fraction of clones closer than r50 / rtot of the cluster
} PhResult;

#define PH_PRE_STEP_MYR 0.5
#define PH_CLUSTER_RV_PRIOR 10.0    // km/s, young clusters follow Galactic rotation closely
#define PH_PULSAR_RV_PRIOR 265.0    // km/s, 1-D neutron-star birth velocity dispersion (Hobbs et al. 2005)
#define PH_STAR_RV_PRIOR 30.0       // km/s, stars without a radial velocity

void ph_init(void);
const Potential *ph_potential(void);

// One draw (rng == NULL: the nominal values). Returns 0, or -1 if no physical draw was found.
int ph_sample_target(const PhTarget *t, McRng *rng, PhaseSpace *out);
int ph_sample_cluster(const PhCluster *c, McRng *rng, PhaseSpace *out);

// Nominal path over tmax_myr plus its spread (from n_spread clones). Free with ph_track_free.
void ph_track_target(const PhTarget *t, double tmax_myr, int n_spread, uint64_t seed, PhTrack *out);
void ph_track_cluster(const PhCluster *c, int n_spread, uint64_t seed, PhTrack *out);
void ph_track_free(PhTrack *t);

// Can the target have come within reach of the cluster? Generous on purpose: margin =
// 3*rtot + 10 pc + 3*(position spread + velocity spread * t).
int ph_prefilter(const PhTrack *target, const PhTrack *cluster, const PhCluster *c);

// The Monte Carlo search for one pair; seed must be a fixed function of the pair (reproducible).
void ph_pair(const PhTarget *t, const PhCluster *c, int n_clones, double dt_myr, uint64_t seed, PhResult *res);

// ---- stage 1 for stars: linear (Gaussian) screening -------------------------------------------
// Stars have small velocity errors, so a draw's deviation from the nominal orbit stays linear in
// its initial deviation: dx(t) = Phi(t) dx0, with Phi the state transition matrix (finite
// differences of 6 neighbouring orbits). The position covariance at time t is then
// A(t) Sigma0 A(t)^T (A = position rows of Phi, Sigma0 = initial 6-D covariance from draws), and a
// pair's relative position is Gaussian with mean m(t) and covariance C(t) = C_target + C_cluster.
// No per-clone integration: ~0.1 s per star against every cluster. Pairs that pass go to ph_pair().
// Not for pulsars: an unknown radial velocity of +-265 km/s is far outside the linear regime.
typedef struct {
    int n;                  // points, PH_PRE_STEP_MYR apart
    Vec3 *pos, *vel;        // nominal, kpc and kpc/Myr
    double (*S)[6];         // position covariance xx xy xz yy yz zz, kpc^2
} PhLin;

typedef struct {
    double chi2;            // min over time of m^T (C + R^2 I)^-1 m, R = ph_screen_r * r50 (or rtot)
    double t;               // Myr ago, where chi2 is smallest
    double d_nom_pc;        // nominal separation then
    double sig_pc;          // RMS position uncertainty of the pair then (per axis)
    double p_r50, p_rtot;   // probability of being inside r50 / rtot at that time
    double vrel;            // nominal relative speed then, km/s
} PhScreen;

void ph_lin_target(const PhTarget *t, double tmax_myr, int n_samp, uint64_t seed, PhLin *out);
// The host half of ph_lin_target (the GPU app integrates on the device): nominal phase space
// (pos.xyz, vel.xyz) and the 6-D covariance of the draws.
void ph_lin_target_init(const PhTarget *t, int n_samp, uint64_t seed, double ps0[6], double sig[6][6]);
#define PH_LIN_EPS_POS 1e-5   // kpc: neighbour-orbit offsets of the state transition matrix
#define PH_LIN_EPS_VEL 1e-5   // kpc/Myr
#define PH_LIN_SUB 10         // leapfrog substeps per PH_PRE_STEP_MYR
void ph_lin_cluster(const PhCluster *c, int n_samp, uint64_t seed, PhLin *out);
void ph_lin_free(PhLin *l);
extern double ph_screen_r;   // screening kernel radius in units of r50 (default 2); 0 = use rtot
void ph_screen(const PhLin *t, const PhLin *c, const PhCluster *cl, uint64_t seed, PhScreen *res);
// ph_screen in two parts: the chi2 scan over time (returns chi2, and where it is smallest), then the
// rest of PhScreen (512-draw probabilities). Callers that drop pairs above a chi2 cut skip the second.
double ph_screen_chi2(const PhLin *t, const PhLin *c, const PhCluster *cl, int *kb, double *sb);
void ph_screen_finish(const PhLin *t, const PhLin *c, const PhCluster *cl, uint64_t seed,
                      double chi2, int kb, double sb, PhScreen *res);

uint64_t ph_seed(const char *target_id, const char *cluster_name);

// File readers (formats: scripts/prepare_clusters.py, scripts/prepare_pulsars.py); return count, -1 on error.
int ph_read_clusters(const char *path, PhCluster **out);
int ph_read_targets(const char *path, PhTarget **out);

#endif
