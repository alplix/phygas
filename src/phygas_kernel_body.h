// Shared device code for the Phygas GPU backends: included by phygas_cuda.cu (as CUDA) and embedded
// as source text by phygas_opencl.c (as OpenCL C). Written in the common subset of both languages:
// PG_GLOBAL / PG_DEVICE / PG_KERNEL / PG_TID / PG_U64 / PG_LONG (64-bit) / PG_INF are defined by the
// includer.
//
// Literal ports of the CPU code (src/phygas.c, src/potential.c, src/montecarlo.[ch]), operation by
// operation, no algebraic rewrites, FP contraction off:
//   pg_lin_k    -- lin_integrate(): nominal orbit + 6 neighbour orbits, position covariance per point
//   pg_chi2_k   -- ph_screen_chi2(): smallest chi2 over time of one (star, cluster) pair
//   pg_finish_k -- ph_screen_finish(): separation, speed and the 512-draw probabilities of one pair
//   pg_pair_k   -- ph_pair_clone(): one Monte Carlo clone of a (target, cluster) pair traced back
//
// Layouts. Star tracks (pg_lin_k output): [star][point][12] = pos.xyz, vel.xyz, S (xx xy xz yy yz zz).
// Cluster tracks: [point][12][cluster] (structure of arrays, so neighbouring threads -- neighbouring
// clusters of one star -- read neighbouring addresses); points past a cluster's own count are unused.

PG_DEVICE void pg_accel(const double *p, double x, double y, double z, double *ax, double *ay, double *az) {
    double R2 = x * x + y * y;
    double R = sqrt(R2);
    double r = sqrt(x * x + y * y + z * z);
    double tx = 0.0, ty = 0.0, tz = 0.0;
    {   // Plummer bulge
        double denom = pow(r * r + p[1] * p[1], 1.5);
        double f = -p[0] / (denom > 0 ? denom : 1e-12);
        tx = tx + x * f; ty = ty + y * f; tz = tz + z * f;
    }
    {   // Miyamoto-Nagai disk
        double zeta = sqrt(z * z + p[4] * p[4]);
        double azeta = p[3] + zeta;
        double D2 = R2 + azeta * azeta;
        double D = sqrt(D2);
        double D3 = D2 * D;
        double aR = -p[2] * R / D3;
        double azz = -p[2] * z * azeta / (zeta > 0 ? zeta : 1e-12) / D3;
        if (R > 1e-12) {
            tx += aR * x / R;
            ty += aR * y / R;
        }
        tz += azz;
    }
    {   // truncated Allen & Santillan halo
        double r_safe = (r > 1e-12 ? r : 1e-12);
        double r_eff = (r_safe < p[7]) ? r_safe : p[7];
        double xx = r_eff / p[6];
        double gm = p[5] * pow(xx, p[8]) / (1.0 + pow(xx, p[8] - 1.0));
        double f = -gm / (r_safe * r_safe * r_safe);
        tx = tx + x * f; ty = ty + y * f; tz = tz + z * f;
    }
    *ax = tx; *ay = ty; *az = tz;
}

// body_step() of src/phygas.c: kick-drift-kick with the acceleration carried over (one force
// evaluation per step). b = pos.xyz, vel.xyz, acc.xyz.
PG_DEVICE void pg_body_step(const double *p, double *b, double dt) {
    double h = 0.5 * dt;
    b[3] = b[3] + b[6] * h; b[4] = b[4] + b[7] * h; b[5] = b[5] + b[8] * h;
    b[0] = b[0] + b[3] * dt; b[1] = b[1] + b[4] * dt; b[2] = b[2] + b[5] * dt;
    pg_accel(p, b[0], b[1], b[2], &b[6], &b[7], &b[8]);
    b[3] = b[3] + b[6] * h; b[4] = b[4] + b[7] * h; b[5] = b[5] + b[8] * h;
}

// Points [k0, k0 + nk) of every star's track. state: 7 bodies x 9 doubles per star, kept between
// launches. ps0: nominal pos.xyz vel.xyz; sig: 6x6 initial covariance (row-major).
PG_KERNEL void pg_lin_k(PG_GLOBAL const double *ps0, PG_GLOBAL const double *sig, PG_GLOBAL double *state,
                      PG_GLOBAL double *trk, PG_GLOBAL const double *params, double eps_pos, double eps_vel,
                      double dt, int sub, int k0, int nk, int npts, int ns) {
    int s = PG_TID;
    if (s >= ns) return;
    double p[9];
    for (int i = 0; i < 9; i++) p[i] = params[i];
    double bd[63];
    if (k0 == 0) {
        for (int j = 0; j < 7; j++) {
            for (int c = 0; c < 6; c++) bd[j * 9 + c] = ps0[s * 6 + c];
            if (j > 0) bd[j * 9 + (j - 1)] += (j - 1) < 3 ? eps_pos : eps_vel;
            pg_accel(p, bd[j * 9], bd[j * 9 + 1], bd[j * 9 + 2], &bd[j * 9 + 6], &bd[j * 9 + 7], &bd[j * 9 + 8]);
        }
    } else {
        for (int i = 0; i < 63; i++) bd[i] = state[(PG_LONG)s * 63 + i];
    }
    double sg[36];
    for (int i = 0; i < 36; i++) sg[i] = sig[(PG_LONG)s * 36 + i];
    for (int k = k0; k < k0 + nk && k < npts; k++) {
        if (k > 0) for (int i = 0; i < sub; i++) for (int j = 0; j < 7; j++) pg_body_step(p, bd + j * 9, dt);
        PG_GLOBAL double *o = trk + ((PG_LONG)s * npts + k) * 12;
        for (int c = 0; c < 6; c++) o[c] = bd[c];
        double A[18];   // 3 x 6
        for (int j = 0; j < 6; j++) {
            double e = j < 3 ? eps_pos : eps_vel;
            A[j] = (bd[(j + 1) * 9] - bd[0]) / e;
            A[6 + j] = (bd[(j + 1) * 9 + 1] - bd[1]) / e;
            A[12 + j] = (bd[(j + 1) * 9 + 2] - bd[2]) / e;
        }
        double AS[18];
        for (int r = 0; r < 3; r++) for (int b = 0; b < 6; b++) {
            double acc = 0;
            for (int a = 0; a < 6; a++) acc += A[r * 6 + a] * sg[a * 6 + b];
            AS[r * 6 + b] = acc;
        }
        double S[9];
        for (int r = 0; r < 3; r++) for (int q = 0; q < 3; q++) {
            double acc = 0;
            for (int b = 0; b < 6; b++) acc += AS[r * 6 + b] * A[q * 6 + b];
            S[r * 3 + q] = acc;
        }
        o[6] = S[0]; o[7] = S[1]; o[8] = S[2]; o[9] = S[4]; o[10] = S[5]; o[11] = S[8];
    }
    for (int i = 0; i < 63; i++) state[(PG_LONG)s * 63 + i] = bd[i];
}

// quad_inv() of src/phygas.c: m^T M^-1 m for symmetric M (xx xy xz yy yz zz); -1 if singular
PG_DEVICE double pg_quad_inv(const double *M, double mx, double my, double mz) {
    double a = M[0], b = M[1], c = M[2], d = M[3], e = M[4], f = M[5];
    double i00 = d * f - e * e, i01 = c * e - b * f, i02 = b * e - c * d;
    double det = a * i00 + b * i01 + c * i02;
    if (!(det > 0)) return -1;
    double i11 = a * f - c * c, i12 = b * c - a * e, i22 = a * d - b * b;
    return (mx * (i00 * mx + i01 * my + i02 * mz) + my * (i01 * mx + i11 * my + i12 * mz)
            + mz * (i02 * mx + i12 * my + i22 * mz)) / det;
}

#define PG_C(k, c, ci) csoa[((PG_LONG)(k) * 12 + (c)) * nc + (ci)]

// Pairs [p0, p0 + np) of the batch (pair = star * nc + cluster; nb * nc < 2^31). Outputs indexed by pair.
PG_KERNEL void pg_chi2_k(PG_GLOBAL const double *trk, PG_GLOBAL const double *csoa, PG_GLOBAL const int *cn,
                       PG_GLOBAL const double *R2s, int nb, int nc, int npts, int p0, int np,
                       PG_GLOBAL double *o_chi2, PG_GLOBAL int *o_kb, PG_GLOBAL double *o_sb) {
    int id = p0 + PG_TID;
    if (PG_TID >= np || id >= nb * nc) return;
    int s = (int)(id / nc), ci = (int)(id % nc);
    int n = npts < cn[ci] ? npts : cn[ci];
    double R2 = R2s[ci];
    PG_GLOBAL const double *t = trk + (PG_LONG)s * npts * 12;
    double best = PG_INF;
    int kb = 0;
    double sb = 0;
    for (int k = 0; k < n; k++) {
        double M[6];
        for (int i = 0; i < 6; i++) M[i] = t[k * 12 + 6 + i] + PG_C(k, 6 + i, ci);
        M[0] += R2; M[3] += R2; M[5] += R2;
        double m0x = t[k * 12] - PG_C(k, 0, ci), m0y = t[k * 12 + 1] - PG_C(k, 1, ci), m0z = t[k * 12 + 2] - PG_C(k, 2, ci);
        double q = pg_quad_inv(M, m0x, m0y, m0z), sv = 0;
        if (k + 1 < n) {
            double m1x = t[(k + 1) * 12] - PG_C(k + 1, 0, ci), m1y = t[(k + 1) * 12 + 1] - PG_C(k + 1, 1, ci),
                   m1z = t[(k + 1) * 12 + 2] - PG_C(k + 1, 2, ci);
            // seg_min(m0, m1, &sv)
            double drx = m1x - m0x, dry = m1y - m0y, drz = m1z - m0z;
            double dd = drx * drx + dry * dry + drz * drz;
            sv = 0;
            if (dd > 0) {
                sv = -(m0x * drx + m0y * dry + m0z * drz) / dd;
                if (sv < 0) sv = 0; else if (sv > 1) sv = 1;
            }
            double msx = m0x + (m1x - m0x) * sv, msy = m0y + (m1y - m0y) * sv, msz = m0z + (m1z - m0z) * sv;
            double qs = pg_quad_inv(M, msx, msy, msz);
            if (qs >= 0 && qs < q) q = qs; else sv = 0;
        }
        if (q >= 0 && q < best) { best = q; kb = k; sb = sv; }
    }
    o_chi2[id] = best;
    o_kb[id] = kb;
    o_sb[id] = sb;
}

// src/montecarlo.h / montecarlo.c
PG_DEVICE PG_U64 pg_rng_next(PG_U64 *st) {
    PG_U64 z = (*st += (PG_U64)0x9E3779B97F4A7C15UL);
    z = (z ^ (z >> 30)) * (PG_U64)0xBF58476D1CE4E5B9UL;
    z = (z ^ (z >> 27)) * (PG_U64)0x94D049BB133111EBUL;
    return z ^ (z >> 31);
}
PG_DEVICE double pg_uniform(PG_U64 *st) { return (double)(pg_rng_next(st) >> 11) * (1.0 / 9007199254740992.0); }
PG_DEVICE double pg_gaussian(PG_U64 *st) {
    double u1 = pg_uniform(st);
    if (u1 <= 0.0) u1 = 1.0 / 18446744073709551616.0;
    double u2 = pg_uniform(st);
    double r = sqrt(-2.0 * log(u1));
    return r * cos(2.0 * 3.14159265358979323846 * u2);
}

// One selected pair per thread, [h0, h0 + nh). out[h * 6 ...]: t, d_nom_pc, sig_pc, p_r50, p_rtot, vrel.
PG_KERNEL void pg_finish_k(PG_GLOBAL const double *trk, PG_GLOBAL const double *csoa, PG_GLOBAL const int *cn,
                         PG_GLOBAL const double *r50s, PG_GLOBAL const double *rts, int nc, int npts,
                         PG_GLOBAL const int *hs, PG_GLOBAL const int *hc, PG_GLOBAL const int *hkb,
                         PG_GLOBAL const double *hsb, PG_GLOBAL const PG_U64 *hseed, int h0, int nh,
                         PG_GLOBAL double *out) {
    int h = h0 + PG_TID;
    if (PG_TID >= nh) return;
    int s = hs[h], ci = hc[h], kb = hkb[h];
    double sb = hsb[h];
    int n = npts < cn[ci] ? npts : cn[ci];
    PG_GLOBAL const double *t = trk + (PG_LONG)s * npts * 12;
    double tm = (kb + sb) * 0.5;   // PH_PRE_STEP_MYR
    int k1 = kb + 1 < n ? kb + 1 : kb;
    double ax = t[kb * 12] - PG_C(kb, 0, ci), ay = t[kb * 12 + 1] - PG_C(kb, 1, ci), az = t[kb * 12 + 2] - PG_C(kb, 2, ci);
    double bx = t[k1 * 12] - PG_C(k1, 0, ci), by = t[k1 * 12 + 1] - PG_C(k1, 1, ci), bz = t[k1 * 12 + 2] - PG_C(k1, 2, ci);
    double mx = ax + (bx - ax) * sb, my = ay + (by - ay) * sb, mz = az + (bz - az) * sb;
    double d_nom = sqrt(mx * mx + my * my + mz * mz) * 1000.0;
    double vx = t[kb * 12 + 3] - PG_C(kb, 3, ci), vy = t[kb * 12 + 4] - PG_C(kb, 4, ci), vz = t[kb * 12 + 5] - PG_C(kb, 5, ci);
    double vrel = sqrt(vx * vx + vy * vy + vz * vz) * 977.7922216807892;
    double C[6];
    for (int i = 0; i < 6; i++) C[i] = t[kb * 12 + 6 + i] + PG_C(kb, 6 + i, ci);
    double sig = sqrt((C[0] + C[3] + C[5]) / 3.0) * 1000.0;
    double l00 = sqrt(fmax(C[0], 0.0)), l10 = l00 > 0 ? C[1] / l00 : 0, l20 = l00 > 0 ? C[2] / l00 : 0;
    double l11 = sqrt(fmax(C[3] - l10 * l10, 0.0)), l21 = l11 > 0 ? (C[4] - l20 * l10) / l11 : 0;
    double l22 = sqrt(fmax(C[5] - l20 * l20 - l21 * l21, 0.0));
    PG_U64 st = hseed[h];
    double r50 = r50s[ci], rt = rts[ci];
    int in50 = 0, intot = 0;
    for (int i = 0; i < 512; i++) {
        double z0 = pg_gaussian(&st), z1 = pg_gaussian(&st), z2 = pg_gaussian(&st);
        double xx = mx + l00 * z0, xy = my + l10 * z0 + l11 * z1, xz = mz + l20 * z0 + l21 * z1 + l22 * z2;
        double r = sqrt(xx * xx + xy * xy + xz * xz);
        if (r < r50) in50++;
        if (r < rt) intot++;
    }
    PG_GLOBAL double *o = out + (PG_LONG)h * 6;
    o[0] = tm; o[1] = d_nom; o[2] = sig;
    o[3] = (double)in50 / 512; o[4] = (double)intot / 512; o[5] = vrel;
}

// ph_pair_clone() of src/phygas.c for clones [0, n): steps (s0, s0 + nk] of each (a clone stops at its
// own step count). ps: 12 per clone (target pos/vel, cluster pos/vel); st: 21 per clone (two bodies
// x pos/vel/acc, then r0); res: 3 per clone (closest distance kpc, its time Myr, speed kpc/Myr).
PG_KERNEL void pg_pair_k(PG_GLOBAL const double *ps, PG_GLOBAL const int *steps, PG_GLOBAL const double *dts,
                         PG_GLOBAL double *st, PG_GLOBAL double *res, PG_GLOBAL const double *params,
                         int s0, int nk, int n) {
    int k = PG_TID;
    if (k >= n) return;
    double p[9];
    for (int i = 0; i < 9; i++) p[i] = params[i];
    double a[9], b[9], r0x, r0y, r0z, best, tbest, vbest;
    PG_GLOBAL double *S = st + (PG_LONG)k * 21;
    PG_GLOBAL double *R = res + (PG_LONG)k * 3;
    if (s0 == 0) {
        for (int c = 0; c < 6; c++) { a[c] = ps[(PG_LONG)k * 12 + c]; b[c] = ps[(PG_LONG)k * 12 + 6 + c]; }
        pg_accel(p, a[0], a[1], a[2], &a[6], &a[7], &a[8]);
        pg_accel(p, b[0], b[1], b[2], &b[6], &b[7], &b[8]);
        r0x = a[0] - b[0]; r0y = a[1] - b[1]; r0z = a[2] - b[2];
        best = sqrt(r0x * r0x + r0y * r0y + r0z * r0z);
        tbest = 0;
        double ux = a[3] - b[3], uy = a[4] - b[4], uz = a[5] - b[5];
        vbest = sqrt(ux * ux + uy * uy + uz * uz);
    } else {
        for (int c = 0; c < 9; c++) { a[c] = S[c]; b[c] = S[9 + c]; }
        r0x = S[18]; r0y = S[19]; r0z = S[20];
        best = R[0]; tbest = R[1]; vbest = R[2];
    }
    int ns = steps[k];
    double dt = dts[k];
    for (int i = s0 + 1; i <= s0 + nk && i <= ns; i++) {
        pg_body_step(p, a, dt);
        pg_body_step(p, b, dt);
        double r1x = a[0] - b[0], r1y = a[1] - b[1], r1z = a[2] - b[2];
        // seg_min(r0, r1, &s)
        double drx = r1x - r0x, dry = r1y - r0y, drz = r1z - r0z;
        double dd = drx * drx + dry * dry + drz * drz, s = 0;
        if (dd > 0) {
            s = -(r0x * drx + r0y * dry + r0z * drz) / dd;
            if (s < 0) s = 0; else if (s > 1) s = 1;
        }
        double mx = r0x + drx * s, my = r0y + dry * s, mz = r0z + drz * s;
        double dm = sqrt(mx * mx + my * my + mz * mz);
        if (dm < best) {
            best = dm; tbest = (i - 1 + s) * -dt;
            double ux = a[3] - b[3], uy = a[4] - b[4], uz = a[5] - b[5];
            vbest = sqrt(ux * ux + uy * uy + uz * uz);
        }
        r0x = r1x; r0y = r1y; r0z = r1z;
    }
    for (int c = 0; c < 9; c++) { S[c] = a[c]; S[9 + c] = b[c]; }
    S[18] = r0x; S[19] = r0y; S[20] = r0z;
    R[0] = best; R[1] = tbest; R[2] = vbest;
}
