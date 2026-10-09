#include "phygas.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "units.h"

static Potential g_pot;
static Vec3 g_sun_pos, g_sun_vel;

void ph_init(void) {
    g_pot = milkyway_potential(&milkyway_potential_default);
    Vec3 zero = {0, 0, 0};
    g_sun_pos = helio_equatorial_to_galactocentric_pos(zero, &solar_params_default);
    g_sun_vel = helio_equatorial_to_galactocentric_vel(zero, &solar_params_default);
}

const Potential *ph_potential(void) { return &g_pot; }

static double gauss(McRng *rng) { return rng ? mc_gaussian(rng) : 0.0; }

// Radial velocity (km/s) a body at this position would have on a circular orbit in the Galactic
// rotation: the prior for clusters and stars without a measured one, and the centre of the pulsar prior.
static double rv_expected(Vec3 pos) {
    double R = hypot(pos.x, pos.y);
    if (R <= 0) return 0.0;
    Vec3 a = potential_acceleration(&g_pot, v3_make(pos.x, pos.y, 0.0));
    double ar = -(a.x * pos.x + a.y * pos.y) / R;
    double vc = ar > 0 ? sqrt(R * ar) : 0.0;
    Vec3 v = v3_make(vc * pos.y / R, -vc * pos.x / R, 0.0);   // the Sun at (-R0,0) moves along +y
    Vec3 los = v3_sub(pos, g_sun_pos);
    double n = v3_norm(los);
    if (n <= 0) return 0.0;
    return nostos_kpcmyr_to_kms(v3_dot(v3_sub(v, g_sun_vel), v3_scale(los, 1.0 / n)));
}

static void finish(Astrometry a, int has_rv, double rv, double e_rv, double prior, McRng *rng, PhaseSpace *out) {
    if (has_rv) {
        a.rv_kms = rv + e_rv * gauss(rng);
    } else {
        a.rv_kms = 0.0;
        PhaseSpace p;
        astrometry_to_galactocentric(a, &solar_params_default, &p);
        a.rv_kms = rv_expected(p.pos) + prior * gauss(rng);
    }
    astrometry_to_galactocentric(a, &solar_params_default, out);
}

int ph_sample_target(const PhTarget *t, McRng *rng, PhaseSpace *out) {
    Astrometry a;
    a.ra_deg = t->ra;
    a.dec_deg = t->dec;
    for (int tries = 0; tries < 50; tries++) {
        double plx = t->plx > 0 ? t->plx + t->e_plx * gauss(rng)
                                : 1000.0 / (t->dist_pc * exp(t->e_lndist * gauss(rng)));
        if (!(plx > 0)) continue;
        a.parallax_mas = plx;
        a.pmra_masyr = t->pmra + t->e_pmra * gauss(rng);
        a.pmdec_masyr = t->pmdec + t->e_pmdec * gauss(rng);
        double prior = t->kind == 'P' ? PH_PULSAR_RV_PRIOR : PH_STAR_RV_PRIOR;
        finish(a, t->has_rv, t->rv, t->e_rv, prior, rng, out);
        return 0;
    }
    return -1;
}

int ph_sample_cluster(const PhCluster *c, McRng *rng, PhaseSpace *out) {
    Astrometry a;
    a.ra_deg = c->ra;
    a.dec_deg = c->dec;
    for (int tries = 0; tries < 50; tries++) {
        double d = c->dist_pc + c->e_dist_pc * gauss(rng);
        if (!(d > 0)) continue;
        a.parallax_mas = 1000.0 / d;
        a.pmra_masyr = c->pmra + c->e_pmra * gauss(rng);
        a.pmdec_masyr = c->pmdec + c->e_pmdec * gauss(rng);
        finish(a, c->has_rv, c->rv, c->e_rv, PH_CLUSTER_RV_PRIOR, rng, out);
        return 0;
    }
    return -1;
}

// ---- integration: kick-drift-kick leapfrog that keeps the acceleration (one evaluation per step)
typedef struct { Vec3 pos, vel, acc; } Body;

static void body_init(Body *b, PhaseSpace ps) {
    b->pos = ps.pos; b->vel = ps.vel; b->acc = potential_acceleration(&g_pot, ps.pos);
}

static inline void body_step(Body *b, double dt) {
    b->vel = v3_add(b->vel, v3_scale(b->acc, 0.5 * dt));
    b->pos = v3_add(b->pos, v3_scale(b->vel, dt));
    b->acc = potential_acceleration(&g_pot, b->pos);
    b->vel = v3_add(b->vel, v3_scale(b->acc, 0.5 * dt));
}

static void track(PhaseSpace ps, double tmax, PhTrack *out) {
    int n = (int)ceil(tmax / PH_PRE_STEP_MYR) + 1;
    out->n = n;
    out->pos = malloc(sizeof(Vec3) * n);
    Body b; body_init(&b, ps);
    const int sub = 10;
    double dt = -PH_PRE_STEP_MYR / sub;
    out->pos[0] = b.pos;
    for (int k = 1; k < n; k++) {
        for (int i = 0; i < sub; i++) body_step(&b, dt);
        out->pos[k] = b.pos;
    }
}

static void spread(PhaseSpace nom, PhaseSpace *draws, int n, PhTrack *out) {
    double sp = 0, sv = 0;
    for (int i = 0; i < n; i++) {
        Vec3 dp = v3_sub(draws[i].pos, nom.pos), dv = v3_sub(draws[i].vel, nom.vel);
        sp += v3_dot(dp, dp); sv += v3_dot(dv, dv);
    }
    out->sig_pos = n ? sqrt(sp / n) : 0;
    out->sig_vel = n ? sqrt(sv / n) : 0;
}

void ph_track_target(const PhTarget *t, double tmax, int n_spread, uint64_t seed, PhTrack *out) {
    PhaseSpace nom, d[64];
    ph_sample_target(t, NULL, &nom);
    McRng rng; mc_rng_seed(&rng, seed);
    if (n_spread > 64) n_spread = 64;
    int k = 0;
    for (int i = 0; i < n_spread; i++) if (ph_sample_target(t, &rng, &d[k]) == 0) k++;
    track(nom, tmax, out);
    spread(nom, d, k, out);
}

void ph_track_cluster(const PhCluster *c, int n_spread, uint64_t seed, PhTrack *out) {
    PhaseSpace nom, d[64];
    ph_sample_cluster(c, NULL, &nom);
    McRng rng; mc_rng_seed(&rng, seed);
    if (n_spread > 64) n_spread = 64;
    int k = 0;
    for (int i = 0; i < n_spread; i++) if (ph_sample_cluster(c, &rng, &d[k]) == 0) k++;
    track(nom, c->tmax_myr, out);
    spread(nom, d, k, out);
}

void ph_track_free(PhTrack *t) { free(t->pos); t->pos = NULL; t->n = 0; }

// smallest |r0 + s (r1 - r0)| for s in [0,1]; *s_out = where
static double seg_min(Vec3 r0, Vec3 r1, double *s_out) {
    Vec3 dr = v3_sub(r1, r0);
    double dd = v3_dot(dr, dr), s = 0;
    if (dd > 0) {
        s = -v3_dot(r0, dr) / dd;
        if (s < 0) s = 0; else if (s > 1) s = 1;
    }
    *s_out = s;
    return v3_norm(v3_add(r0, v3_scale(dr, s)));
}

int ph_prefilter(const PhTrack *tt, const PhTrack *ct, const PhCluster *c) {
    int n = tt->n < ct->n ? tt->n : ct->n;
    double sp = hypot(tt->sig_pos, ct->sig_pos), sv = hypot(tt->sig_vel, ct->sig_vel);
    double r0 = (3.0 * c->rtot_pc + 10.0) / 1000.0;
    for (int k = 0; k + 1 < n; k++) {
        double s, d = seg_min(v3_sub(tt->pos[k], ct->pos[k]), v3_sub(tt->pos[k + 1], ct->pos[k + 1]), &s);
        double t = (k + s) * PH_PRE_STEP_MYR;
        if (d < r0 + 3.0 * (sp + sv * t)) return 1;
    }
    return 0;
}

static int cmp_double(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b;
    return x < y ? -1 : x > y;
}

static double quantile(const double *v, int n, double q) {
    double pos = q * (n - 1);
    int i = (int)pos;
    if (i >= n - 1) return v[n - 1];
    return v[i] + (pos - i) * (v[i + 1] - v[i]);
}

int ph_pair_draw(const PhTarget *t, const PhCluster *c, int n_clones, double dt_myr, uint64_t seed,
                 double *ps, int *steps_out, double *dt_out) {
    McRng rng; mc_rng_seed(&rng, seed);
    int steps = (int)ceil(c->tmax_myr / dt_myr);
    *steps_out = steps;
    *dt_out = -c->tmax_myr / steps;
    int m = 0;
    for (int k = 0; k < n_clones; k++) {
        PhaseSpace ps_t, ps_c;
        if (ph_sample_target(t, &rng, &ps_t) || ph_sample_cluster(c, &rng, &ps_c)) continue;
        double *o = ps + (size_t)m * 12;
        o[0] = ps_t.pos.x; o[1] = ps_t.pos.y; o[2] = ps_t.pos.z; o[3] = ps_t.vel.x; o[4] = ps_t.vel.y; o[5] = ps_t.vel.z;
        o[6] = ps_c.pos.x; o[7] = ps_c.pos.y; o[8] = ps_c.pos.z; o[9] = ps_c.vel.x; o[10] = ps_c.vel.y; o[11] = ps_c.vel.z;
        m++;
    }
    return m;
}

void ph_pair_clone(const double *ps, int steps, double dt, double *best_out, double *tbest_out, double *vbest_out) {
    PhaseSpace pt = {{ps[0], ps[1], ps[2]}, {ps[3], ps[4], ps[5]}}, pc = {{ps[6], ps[7], ps[8]}, {ps[9], ps[10], ps[11]}};
    Body a, b; body_init(&a, pt); body_init(&b, pc);
    Vec3 r0 = v3_sub(a.pos, b.pos);
    double best = v3_norm(r0), tbest = 0, vbest = v3_norm(v3_sub(a.vel, b.vel));
    for (int i = 1; i <= steps; i++) {
        body_step(&a, dt); body_step(&b, dt);
        Vec3 r1 = v3_sub(a.pos, b.pos);
        double s, dm = seg_min(r0, r1, &s);
        if (dm < best) {
            best = dm; tbest = (i - 1 + s) * -dt;
            vbest = v3_norm(v3_sub(a.vel, b.vel));
        }
        r0 = r1;
    }
    *best_out = best; *tbest_out = tbest; *vbest_out = vbest;
}

void ph_pair_reduce(const PhCluster *c, int m, const double *best, const double *tbest, const double *vbest, PhResult *res) {
    double *d = malloc(sizeof(double) * (m > 0 ? m : 1) * 3), *tt = d + (m > 0 ? m : 1), *vv = d + 2 * (m > 0 ? m : 1);
    int in50 = 0, intot = 0;
    for (int k = 0; k < m; k++) {
        d[k] = best[k] * 1000.0; tt[k] = tbest[k]; vv[k] = nostos_kpcmyr_to_kms(vbest[k]);
        if (d[k] < c->r50_pc) in50++;
        if (d[k] < c->rtot_pc) intot++;
    }
    memset(res, 0, sizeof(*res));
    res->n = m;
    if (m > 0) {
        res->f_r50 = (double)in50 / m; res->f_rtot = (double)intot / m;
        qsort(d, m, sizeof(double), cmp_double);
        qsort(tt, m, sizeof(double), cmp_double);
        qsort(vv, m, sizeof(double), cmp_double);
        res->d05 = quantile(d, m, 0.05); res->d50 = quantile(d, m, 0.5); res->d95 = quantile(d, m, 0.95);
        res->t05 = quantile(tt, m, 0.05); res->t50 = quantile(tt, m, 0.5); res->t95 = quantile(tt, m, 0.95);
        res->v05 = quantile(vv, m, 0.05); res->v50 = quantile(vv, m, 0.5); res->v95 = quantile(vv, m, 0.95);
    }
    free(d);
}

void ph_pair(const PhTarget *t, const PhCluster *c, int n_clones, double dt_myr, uint64_t seed, PhResult *res) {
    double *ps = malloc(sizeof(double) * 12 * (n_clones > 0 ? n_clones : 1));
    double *out = malloc(sizeof(double) * 3 * (n_clones > 0 ? n_clones : 1));
    int steps;
    double dt;
    int m = ph_pair_draw(t, c, n_clones, dt_myr, seed, ps, &steps, &dt);
    for (int k = 0; k < m; k++) ph_pair_clone(ps + (size_t)k * 12, steps, dt, &out[k], &out[m + k], &out[2 * m + k]);
    ph_pair_reduce(c, m, out, out + m, out + 2 * m, res);
    free(ps); free(out);
}

// ---- stage 1: linear screening ---------------------------------------------------------------
#define PH_EPS_POS PH_LIN_EPS_POS     // kpc (0.01 pc)
#define PH_EPS_VEL PH_LIN_EPS_VEL     // kpc/Myr (~0.01 km/s)

// initial 6-D covariance from the draws
static void lin_cov(PhaseSpace nom, const PhaseSpace *draws, int nd, double sig[6][6]) {
    memset(sig, 0, sizeof(double[6][6]));
    for (int i = 0; i < nd; i++) {
        double d[6] = {draws[i].pos.x - nom.pos.x, draws[i].pos.y - nom.pos.y, draws[i].pos.z - nom.pos.z,
                       draws[i].vel.x - nom.vel.x, draws[i].vel.y - nom.vel.y, draws[i].vel.z - nom.vel.z};
        for (int a = 0; a < 6; a++) for (int b = 0; b < 6; b++) sig[a][b] += d[a] * d[b];
    }
    for (int a = 0; a < 6; a++) for (int b = 0; b < 6; b++) sig[a][b] /= nd > 0 ? nd : 1;
}

static void lin_integrate(PhaseSpace nom, double sig[6][6], double tmax, PhLin *out);

static void lin_build(PhaseSpace nom, PhaseSpace *draws, int nd, double tmax, PhLin *out) {
    double sig[6][6];
    lin_cov(nom, draws, nd, sig);
    lin_integrate(nom, sig, tmax, out);
}

void ph_lin_target_init(const PhTarget *t, int n_samp, uint64_t seed, double ps0[6], double sig[6][6]) {
    PhaseSpace nom, *d = malloc(sizeof(PhaseSpace) * n_samp);
    ph_sample_target(t, NULL, &nom);
    McRng rng; mc_rng_seed(&rng, seed);
    int k = 0;
    for (int i = 0; i < n_samp; i++) if (ph_sample_target(t, &rng, &d[k]) == 0) k++;
    lin_cov(nom, d, k, sig);
    free(d);
    ps0[0] = nom.pos.x; ps0[1] = nom.pos.y; ps0[2] = nom.pos.z;
    ps0[3] = nom.vel.x; ps0[4] = nom.vel.y; ps0[5] = nom.vel.z;
}

static void lin_integrate(PhaseSpace nom, double sig[6][6], double tmax, PhLin *out) {
    // nominal orbit + 6 neighbours, one per phase-space axis
    Body bd[7];
    body_init(&bd[0], nom);
    for (int j = 0; j < 6; j++) {
        PhaseSpace p = nom;
        double e = j < 3 ? PH_EPS_POS : PH_EPS_VEL;
        if (j == 0) p.pos.x += e; else if (j == 1) p.pos.y += e; else if (j == 2) p.pos.z += e;
        else if (j == 3) p.vel.x += e; else if (j == 4) p.vel.y += e; else p.vel.z += e;
        body_init(&bd[j + 1], p);
    }
    int n = (int)ceil(tmax / PH_PRE_STEP_MYR) + 1;
    out->n = n;
    out->pos = malloc(sizeof(Vec3) * n);
    out->vel = malloc(sizeof(Vec3) * n);
    out->S = malloc(sizeof(double[6]) * n);
    const int sub = PH_LIN_SUB;
    double dt = -PH_PRE_STEP_MYR / sub;
    for (int k = 0; k < n; k++) {
        if (k > 0) for (int i = 0; i < sub; i++) for (int j = 0; j < 7; j++) body_step(&bd[j], dt);
        out->pos[k] = bd[0].pos;
        out->vel[k] = bd[0].vel;
        double A[3][6];
        for (int j = 0; j < 6; j++) {
            double e = j < 3 ? PH_EPS_POS : PH_EPS_VEL;
            Vec3 d = v3_sub(bd[j + 1].pos, bd[0].pos);
            A[0][j] = d.x / e; A[1][j] = d.y / e; A[2][j] = d.z / e;
        }
        double AS[3][6] = {{0}};
        for (int r = 0; r < 3; r++) for (int b = 0; b < 6; b++)
            for (int a = 0; a < 6; a++) AS[r][b] += A[r][a] * sig[a][b];
        double S[3][3] = {{0}};
        for (int r = 0; r < 3; r++) for (int s = 0; s < 3; s++)
            for (int b = 0; b < 6; b++) S[r][s] += AS[r][b] * A[s][b];
        out->S[k][0] = S[0][0]; out->S[k][1] = S[0][1]; out->S[k][2] = S[0][2];
        out->S[k][3] = S[1][1]; out->S[k][4] = S[1][2]; out->S[k][5] = S[2][2];
    }
}

void ph_lin_target(const PhTarget *t, double tmax, int n_samp, uint64_t seed, PhLin *out) {
    PhaseSpace nom, *d = malloc(sizeof(PhaseSpace) * n_samp);
    ph_sample_target(t, NULL, &nom);
    McRng rng; mc_rng_seed(&rng, seed);
    int k = 0;
    for (int i = 0; i < n_samp; i++) if (ph_sample_target(t, &rng, &d[k]) == 0) k++;
    lin_build(nom, d, k, tmax, out);
    free(d);
}

void ph_lin_cluster(const PhCluster *c, int n_samp, uint64_t seed, PhLin *out) {
    PhaseSpace nom, *d = malloc(sizeof(PhaseSpace) * n_samp);
    ph_sample_cluster(c, NULL, &nom);
    McRng rng; mc_rng_seed(&rng, seed);
    int k = 0;
    for (int i = 0; i < n_samp; i++) if (ph_sample_cluster(c, &rng, &d[k]) == 0) k++;
    lin_build(nom, d, k, c->tmax_myr, out);
    free(d);
}

void ph_lin_free(PhLin *l) { free(l->pos); free(l->vel); free(l->S); l->pos = l->vel = NULL; l->S = NULL; l->n = 0; }

// symmetric 3x3 (xx xy xz yy yz zz): inverse quadratic form m^T M^-1 m; returns -1 if singular
static double quad_inv(const double M[6], Vec3 m) {
    double a = M[0], b = M[1], c = M[2], d = M[3], e = M[4], f = M[5];
    double i00 = d * f - e * e, i01 = c * e - b * f, i02 = b * e - c * d;
    double det = a * i00 + b * i01 + c * i02;
    if (!(det > 0)) return -1;
    double i11 = a * f - c * c, i12 = b * c - a * e, i22 = a * d - b * b;
    return (m.x * (i00 * m.x + i01 * m.y + i02 * m.z) + m.y * (i01 * m.x + i11 * m.y + i12 * m.z)
            + m.z * (i02 * m.x + i12 * m.y + i22 * m.z)) / det;
}

double ph_screen_r = 2.0;

double ph_screen_chi2(const PhLin *tl, const PhLin *cl, const PhCluster *c, int *kb_out, double *sb_out) {
    int n = tl->n < cl->n ? tl->n : cl->n;
    // kernel radius for "inside the cluster": PH_SCREEN_R (r50 multiples; 0 = total radius)
    double R = (ph_screen_r > 0 ? ph_screen_r * c->r50_pc : c->rtot_pc) / 1000.0, R2 = R * R;
    double best = INFINITY;
    int kb = 0;
    double sb = 0;
    for (int k = 0; k < n; k++) {
        double M[6];
        for (int i = 0; i < 6; i++) M[i] = tl->S[k][i] + cl->S[k][i];
        M[0] += R2; M[3] += R2; M[5] += R2;
        Vec3 m0 = v3_sub(tl->pos[k], cl->pos[k]);
        double q = quad_inv(M, m0), s = 0;
        if (k + 1 < n) {   // the means move ~linearly within a step: check their closest point too
            Vec3 m1 = v3_sub(tl->pos[k + 1], cl->pos[k + 1]);
            double dmin = seg_min(m0, m1, &s);
            (void)dmin;
            Vec3 ms = v3_add(m0, v3_scale(v3_sub(m1, m0), s));
            double qs = quad_inv(M, ms);
            if (qs >= 0 && qs < q) q = qs; else s = 0;
        }
        if (q >= 0 && q < best) { best = q; kb = k; sb = s; }
    }
    *kb_out = kb; *sb_out = sb;
    return best;
}

void ph_screen_finish(const PhLin *tl, const PhLin *cl, const PhCluster *c, uint64_t seed,
                      double best, int kb, double sb, PhScreen *res) {
    int n = tl->n < cl->n ? tl->n : cl->n;
    memset(res, 0, sizeof(*res));
    res->chi2 = best;
    res->t = (kb + sb) * PH_PRE_STEP_MYR;
    int k1 = kb + 1 < n ? kb + 1 : kb;
    Vec3 m = v3_add(v3_sub(tl->pos[kb], cl->pos[kb]),
                    v3_scale(v3_sub(v3_sub(tl->pos[k1], cl->pos[k1]), v3_sub(tl->pos[kb], cl->pos[kb])), sb));
    res->d_nom_pc = v3_norm(m) * 1000.0;
    res->vrel = nostos_kpcmyr_to_kms(v3_norm(v3_sub(tl->vel[kb], cl->vel[kb])));
    double C[6];
    for (int i = 0; i < 6; i++) C[i] = tl->S[kb][i] + cl->S[kb][i];
    res->sig_pc = sqrt((C[0] + C[3] + C[5]) / 3.0) * 1000.0;
    // probability of being inside r50 / rtot at that time: 512 Gaussian draws (Cholesky of C)
    double l00 = sqrt(fmax(C[0], 0)), l10 = l00 > 0 ? C[1] / l00 : 0, l20 = l00 > 0 ? C[2] / l00 : 0;
    double l11 = sqrt(fmax(C[3] - l10 * l10, 0)), l21 = l11 > 0 ? (C[4] - l20 * l10) / l11 : 0;
    double l22 = sqrt(fmax(C[5] - l20 * l20 - l21 * l21, 0));
    McRng rng; mc_rng_seed(&rng, seed);
    double r50 = c->r50_pc / 1000.0, rt = c->rtot_pc / 1000.0;
    int in50 = 0, intot = 0;
    const int N = 512;
    for (int i = 0; i < N; i++) {
        double z0 = mc_gaussian(&rng), z1 = mc_gaussian(&rng), z2 = mc_gaussian(&rng);
        Vec3 x = v3_make(m.x + l00 * z0, m.y + l10 * z0 + l11 * z1, m.z + l20 * z0 + l21 * z1 + l22 * z2);
        double r = v3_norm(x);
        if (r < r50) in50++;
        if (r < rt) intot++;
    }
    res->p_r50 = (double)in50 / N;
    res->p_rtot = (double)intot / N;
}

void ph_screen(const PhLin *tl, const PhLin *cl, const PhCluster *c, uint64_t seed, PhScreen *res) {
    int kb;
    double sb, best = ph_screen_chi2(tl, cl, c, &kb, &sb);
    ph_screen_finish(tl, cl, c, seed, best, kb, sb, res);
}

uint64_t ph_seed(const char *a, const char *b) {
    uint64_t h = 1469598103934665603ULL;   // FNV-1a over "a|b"
    for (const char *p = a; *p; p++) { h ^= (unsigned char)*p; h *= 1099511628211ULL; }
    h ^= '|'; h *= 1099511628211ULL;
    for (const char *p = b; *p; p++) { h ^= (unsigned char)*p; h *= 1099511628211ULL; }
    return h;
}

static double parse_num(const char *s) { return strcmp(s, "nan") == 0 ? NAN : atof(s); }

int ph_read_clusters(const char *path, PhCluster **out) {
    FILE *f = fopen(path, "r");
    if (!f) return -1;
    int cap = 1024, n = 0;
    PhCluster *v = malloc(sizeof(PhCluster) * cap);
    char line[1024], f_rv[32], f_erv[32], kind[4];
    int idx;
    while (fgets(line, sizeof line, f)) {
        if (n == cap) { cap *= 2; v = realloc(v, sizeof(PhCluster) * cap); }
        PhCluster *c = &v[n];
        if (sscanf(line, "%d %39s %3s %lf %lf %lf %lf %lf %lf %lf %lf %31s %31s %lf %lf %lf %lf %lf %lf",
                   &idx, c->name, kind, &c->ra, &c->dec, &c->dist_pc, &c->e_dist_pc, &c->pmra, &c->e_pmra,
                   &c->pmdec, &c->e_pmdec, f_rv, f_erv, &c->r50_pc, &c->rtot_pc, &c->age16, &c->age50,
                   &c->age84, &c->tmax_myr) != 19) continue;
        c->kind = kind[0];
        c->rv = parse_num(f_rv); c->e_rv = parse_num(f_erv);
        c->has_rv = !isnan(c->rv) && !isnan(c->e_rv);
        n++;
    }
    fclose(f);
    *out = v;
    return n;
}

int ph_read_targets(const char *path, PhTarget **out) {
    FILE *f = fopen(path, "r");
    if (!f) return -1;
    int cap = 1024, n = 0;
    PhTarget *v = malloc(sizeof(PhTarget) * cap);
    char line[1024], f_rv[32], f_erv[32], kind[4];
    while (fgets(line, sizeof line, f)) {
        if (n == cap) { cap *= 2; v = realloc(v, sizeof(PhTarget) * cap); }
        PhTarget *t = &v[n];
        if (sscanf(line, "%39s %3s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %31s %31s",
                   t->id, kind, &t->ra, &t->dec, &t->plx, &t->e_plx, &t->dist_pc, &t->e_lndist,
                   &t->pmra, &t->e_pmra, &t->pmdec, &t->e_pmdec, f_rv, f_erv) != 14) continue;
        t->kind = kind[0];
        t->rv = parse_num(f_rv); t->e_rv = parse_num(f_erv);
        t->has_rv = !isnan(t->rv) && !isnan(t->e_rv);
        n++;
    }
    fclose(f);
    *out = v;
    return n;
}
