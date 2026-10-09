// Phygas BOINC app. Two input files:
//   "clusters" -- the young-cluster list (scripts/prepare_clusters.py), the same for every task
//                 (sticky on the volunteer's machine, downloaded once)
//   "in"       -- line 1: "PHYGAS 1 screen <chi2_max> <n_samp> <r_kernel> <p_min>"
//                      or "PHYGAS 1 pair <n_clones> <dt_myr>"
//                 then one target per line (src/phygas.h PhTarget format). In pair mode a target
//                 line may end with "clusters=NAME1,NAME2,..." to run only those clusters (stage 2
//                 for stars); without it every cluster that passes the prefilter is run (pulsars).
// Output "out": one line per pair, the same columns as phygas_screen / phygas_cli, then
// "# done <n_targets>". Checkpoint: number of finished targets (output is appended).
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "boinc_api.h"
#include "filesys.h"

extern "C" {
#include "phygas.h"
}

static int fail(const char *msg) {
    fprintf(stderr, "%s\n", msg);
    boinc_finish(1);
    return 1;
}

#ifdef PG_GPU
// GPU build (-DPG_GPU with phygas_cuda.cu or phygas_opencl.c): screen mode runs on the GPU
// (phygas_gpu.h); pair mode (pulsars, stage 2) is the CPU code below. The host keeps everything
// that sets the numbers' inputs -- cluster tracks, star draws and covariances, pair seeds -- as the
// CPU app computes them; the device runs lin_integrate / ph_screen_chi2 / ph_screen_finish ported
// operation by operation, so the output matches the CPU app's to rounding.
#include <vector>
#include "phygas_gpu.h"

static int env_int(const char *name, int def) {
    const char *v = getenv(name);
    int x = v ? atoi(v) : 0;
    return x > 0 ? x : def;
}

struct DevBufs {
    PgBuf ps0 = 0, sig = 0, state = 0, trk = 0, chi2 = 0, kb = 0, sb = 0;
    void release() { pg_free(ps0); pg_free(sig); pg_free(state); pg_free(trk); pg_free(chi2); pg_free(kb); pg_free(sb); *this = DevBufs(); }
    bool alloc(int B, int npts, int nc) {
        ps0 = pg_alloc(sizeof(double) * 6 * B); sig = pg_alloc(sizeof(double) * 36 * B);
        state = pg_alloc(sizeof(double) * 63 * B); trk = pg_alloc(sizeof(double) * 12 * (size_t)npts * B);
        chi2 = pg_alloc(sizeof(double) * (size_t)B * nc); kb = pg_alloc(sizeof(int) * (size_t)B * nc);
        sb = pg_alloc(sizeof(double) * (size_t)B * nc);
        if (ps0 && sig && state && trk && chi2 && kb && sb) return true;
        release();
        return false;
    }
};

// Screen targets [start, nt) and append the lines the CPU loop would print. 0 = ok.
static int gpu_screen(const PhTarget *tg, int nt, int start, const PhCluster *cl, int nc, const PhLin *clin,
                      double tmax, double chi2_max, int n_samp, double p_min, FILE *out) {
    const int npts = (int)ceil(tmax / PH_PRE_STEP_MYR) + 1;   // lin_integrate's point count
    // clusters: tracks as structure of arrays, kernel radii
    std::vector<double> csoa((size_t)npts * 12 * nc, 0.0), R2s(nc), r50s(nc), rts(nc);
    std::vector<int> cn(nc);
    for (int ci = 0; ci < nc; ci++) {
        const PhLin &L = clin[ci];
        cn[ci] = L.n;
        for (int k = 0; k < L.n && k < npts; k++) {
            const double v[12] = {L.pos[k].x, L.pos[k].y, L.pos[k].z, L.vel[k].x, L.vel[k].y, L.vel[k].z,
                                  L.S[k][0], L.S[k][1], L.S[k][2], L.S[k][3], L.S[k][4], L.S[k][5]};
            for (int c = 0; c < 12; c++) csoa[((size_t)k * 12 + c) * nc + ci] = v[c];
        }
        double R = (ph_screen_r > 0 ? ph_screen_r * cl[ci].r50_pc : cl[ci].rtot_pc) / 1000.0;
        R2s[ci] = R * R;
        r50s[ci] = cl[ci].r50_pc / 1000.0;
        rts[ci] = cl[ci].rtot_pc / 1000.0;
    }
    const MilkyWayParams &mw = milkyway_potential_default;
    const double prm[9] = {mw.bulge_GM, mw.bulge_b, mw.disk_GM, mw.disk_a, mw.disk_b,
                           mw.halo_GM, mw.halo_a, mw.halo_cutoff, mw.halo_gamma};
    PgBuf d_csoa = pg_alloc(sizeof(double) * csoa.size()), d_cn = pg_alloc(sizeof(int) * nc),
          d_R2 = pg_alloc(sizeof(double) * nc), d_r50 = pg_alloc(sizeof(double) * nc),
          d_rt = pg_alloc(sizeof(double) * nc), d_prm = pg_alloc(sizeof prm);
    if (!d_csoa || !d_cn || !d_R2 || !d_r50 || !d_rt || !d_prm ||
        pg_put(d_csoa, csoa.data(), sizeof(double) * csoa.size()) || pg_put(d_cn, cn.data(), sizeof(int) * nc) ||
        pg_put(d_R2, R2s.data(), sizeof(double) * nc) || pg_put(d_r50, r50s.data(), sizeof(double) * nc) ||
        pg_put(d_rt, rts.data(), sizeof(double) * nc) || pg_put(d_prm, prm, sizeof prm)) {
        fprintf(stderr, "can't set up the cluster tracks on the GPU\n");
        return 1;
    }
    // stars per batch: PHYGAS_BATCH_STARS (default 1024), halved until the buffers fit
    int B = env_int("PHYGAS_BATCH_STARS", 1024);
    DevBufs d;
    while (!d.alloc(B, npts, nc)) {
        if (B <= 16) { fprintf(stderr, "not enough GPU memory\n"); return 1; }
        B /= 2;
    }
    fprintf(stderr, "GPU screening: %d stars per batch, %d clusters, %d points per track\n", B, nc, npts);
    const double dt = -PH_PRE_STEP_MYR / PH_LIN_SUB;
    std::vector<double> ps0((size_t)6 * B), sig((size_t)36 * B), chi2((size_t)B * nc), sb((size_t)B * nc);
    std::vector<int> kb((size_t)B * nc);
    int lin_slice = 4, chi2_slice = 65536, fin_slice = 16384;
    double t_start = pg_now();
    int rc = 0;
    for (int j0 = start; j0 < nt && !rc; j0 += B) {
        const int nb = nt - j0 < B ? nt - j0 : B;
        for (int s = 0; s < nb; s++) {
            double S6[6][6];
            ph_lin_target_init(&tg[j0 + s], n_samp, ph_seed("lin", tg[j0 + s].id), &ps0[6 * s], S6);
            memcpy(&sig[36 * s], S6, sizeof S6);
        }
        if (pg_put(d.ps0, ps0.data(), sizeof(double) * 6 * nb) || pg_put(d.sig, sig.data(), sizeof(double) * 36 * nb)) { rc = 1; break; }
        for (int k = 0; k < npts && !rc; ) {
            int nk = npts - k < lin_slice ? npts - k : lin_slice;
            double t0 = pg_now();
            rc = pg_lin(d.ps0, d.sig, d.state, d.trk, d_prm, PH_LIN_EPS_POS, PH_LIN_EPS_VEL, dt, PH_LIN_SUB, k, nk, npts, nb);
            k += nk;
            lin_slice = pg_next_slice(nk, pg_now() - t0, npts);
        }
        const int total = nb * nc;
        for (int p = 0; p < total && !rc; ) {
            int np = total - p < chi2_slice ? total - p : chi2_slice;
            double t0 = pg_now();
            rc = pg_chi2(d.trk, d_csoa, d_cn, d_R2, nb, nc, npts, p, np, d.chi2, d.kb, d.sb);
            p += np;
            chi2_slice = pg_next_slice(np, pg_now() - t0, total);
        }
        if (rc || pg_get(d.chi2, chi2.data(), sizeof(double) * total) || pg_get(d.kb, kb.data(), sizeof(int) * total) ||
            pg_get(d.sb, sb.data(), sizeof(double) * total)) { rc = 1; break; }
        // the pairs the CPU loop keeps past the chi2 cut, in its order (star, then cluster)
        std::vector<int> hs, hc, hkb;
        std::vector<double> hsb, hchi2;
        std::vector<uint64_t> hseed;
        for (int s = 0; s < nb; s++) for (int ci = 0; ci < nc; ci++) {
            size_t i = (size_t)s * nc + ci;
            if (!(chi2[i] < chi2_max)) continue;
            hs.push_back(s); hc.push_back(ci); hkb.push_back(kb[i]); hsb.push_back(sb[i]); hchi2.push_back(chi2[i]);
            hseed.push_back(ph_seed(tg[j0 + s].id, cl[ci].name));
        }
        const int nh = (int)hs.size();
        std::vector<double> res((size_t)6 * nh);
        if (nh > 0) {
            PgBuf b_hs = pg_alloc(sizeof(int) * nh), b_hc = pg_alloc(sizeof(int) * nh), b_hkb = pg_alloc(sizeof(int) * nh),
                  b_hsb = pg_alloc(sizeof(double) * nh), b_seed = pg_alloc(sizeof(uint64_t) * nh),
                  b_out = pg_alloc(sizeof(double) * 6 * nh);
            rc = !b_hs || !b_hc || !b_hkb || !b_hsb || !b_seed || !b_out ||
                 pg_put(b_hs, hs.data(), sizeof(int) * nh) || pg_put(b_hc, hc.data(), sizeof(int) * nh) ||
                 pg_put(b_hkb, hkb.data(), sizeof(int) * nh) || pg_put(b_hsb, hsb.data(), sizeof(double) * nh) ||
                 pg_put(b_seed, hseed.data(), sizeof(uint64_t) * nh);
            for (int h = 0; h < nh && !rc; ) {
                int n = nh - h < fin_slice ? nh - h : fin_slice;
                double t0 = pg_now();
                rc = pg_finish(d.trk, d_csoa, d_cn, d_r50, d_rt, nc, npts, b_hs, b_hc, b_hkb, b_hsb, b_seed, h, n, b_out);
                h += n;
                fin_slice = pg_next_slice(n, pg_now() - t0, nh);
            }
            if (!rc) rc = pg_get(b_out, res.data(), sizeof(double) * 6 * nh);
            pg_free(b_hs); pg_free(b_hc); pg_free(b_hkb); pg_free(b_hsb); pg_free(b_seed); pg_free(b_out);
            if (rc) break;
        }
        for (int h = 0; h < nh; h++) {
            const double *r = &res[(size_t)6 * h];   // t, d_nom_pc, sig_pc, p_r50, p_rtot, vrel
            if (r[3] < p_min) continue;
            const PhCluster &c = cl[hc[h]];
            fprintf(out, "%s %s %c %.3f %.3f %.2f %.2f %.4f %.4f %.2f\n", tg[j0 + hs[h]].id, c.name, c.kind,
                    hchi2[h], r[0], r[1], r[2], r[3], r[4], r[5]);
        }
        fflush(out);
        FILE *cw = boinc_fopen("phygas_checkpoint", "w");
        if (cw) { fprintf(cw, "%d\n", j0 + nb); fclose(cw); }
        boinc_checkpoint_completed();
        boinc_fraction_done((double)(j0 + nb) / nt);
    }
    if (!rc) fprintf(stderr, "GPU screening: %d stars in %.1f s\n", nt - start, pg_now() - t_start);
    d.release();
    pg_free(d_csoa); pg_free(d_cn); pg_free(d_R2); pg_free(d_r50); pg_free(d_rt); pg_free(d_prm);
    return rc;
}

static bool in_list(const char *name, const char *list);

// Pair mode (pulsars; stage 2 for stars): the clones of every pair are drawn on the host exactly as
// ph_pair draws them (ph_pair_draw), traced back on the device (pg_pair_k = ph_pair_clone) and
// reduced on the host (ph_pair_reduce), so every line equals the CPU app's to rounding. Clones of many
// pairs and targets go to the GPU together (PHYGAS_BATCH_CLONES, default 262144); output and the
// checkpoint advance only by whole targets.
static int gpu_pair(const PhTarget *tg, char *const *only, int nt, int start, const PhCluster *cl, int nc,
                    const PhTrack *ctr, double tmax, int n_clones, double dt_myr, FILE *out) {
    const MilkyWayParams &mw = milkyway_potential_default;
    const double prm[9] = {mw.bulge_GM, mw.bulge_b, mw.disk_GM, mw.disk_a, mw.disk_b,
                           mw.halo_GM, mw.halo_a, mw.halo_cutoff, mw.halo_gamma};
    int batch = env_int("PHYGAS_BATCH_CLONES", 262144);
    if (batch < n_clones) batch = n_clones;
    PgBuf d_prm = pg_alloc(sizeof prm), d_ps = 0, d_steps = 0, d_dts = 0, d_st = 0, d_res = 0;
    int cap = 0;
    for (cap = 4 * batch; cap >= n_clones; cap /= 2) {   // device room for up to cap clones
        d_ps = pg_alloc(sizeof(double) * 12 * (size_t)cap); d_steps = pg_alloc(sizeof(int) * (size_t)cap);
        d_dts = pg_alloc(sizeof(double) * (size_t)cap); d_st = pg_alloc(sizeof(double) * 21 * (size_t)cap);
        d_res = pg_alloc(sizeof(double) * 3 * (size_t)cap);
        if (d_ps && d_steps && d_dts && d_st && d_res) break;
        pg_free(d_ps); pg_free(d_steps); pg_free(d_dts); pg_free(d_st); pg_free(d_res);
        d_ps = d_steps = d_dts = d_st = d_res = 0;
    }
    if (!d_prm || !d_ps || pg_put(d_prm, prm, sizeof prm)) { fprintf(stderr, "not enough GPU memory\n"); return 1; }
    fprintf(stderr, "GPU pairs: %d clones per pair, up to %d clones per batch\n", n_clones, cap);

    struct Pend { int j, ci, m; size_t off; };
    struct Done { int j, ci; PhResult r; };
    std::vector<Pend> pend;
    std::vector<Done> done;
    std::vector<double> cps, cdts, res;
    std::vector<int> csteps;
    cps.reserve((size_t)12 * cap);
    int slice = 20, rc = 0, written = start;
    long long n_pairs = 0, n_cl = 0;
    double t_start = pg_now();

    auto flush = [&]() -> int {
        const int n = (int)csteps.size();
        if (n > 0) {
            int smax = 0;
            for (int s : csteps) if (s > smax) smax = s;
            if (pg_put(d_ps, cps.data(), sizeof(double) * 12 * n) || pg_put(d_steps, csteps.data(), sizeof(int) * n) ||
                pg_put(d_dts, cdts.data(), sizeof(double) * n)) return 1;
            for (int s0 = 0; s0 < smax; ) {
                int k = smax - s0 < slice ? smax - s0 : slice;
                double t0 = pg_now();
                if (pg_pair(d_ps, d_steps, d_dts, d_st, d_res, d_prm, s0, k, n)) return 1;
                s0 += k;
                slice = pg_next_slice(k, pg_now() - t0, smax);
            }
            res.resize((size_t)3 * n);
            if (pg_get(d_res, res.data(), sizeof(double) * 3 * n)) return 1;
        }
        std::vector<double> b, t, v;
        for (const Pend &p : pend) {
            b.resize(p.m > 0 ? p.m : 1); t.resize(b.size()); v.resize(b.size());
            for (int k = 0; k < p.m; k++) {
                const double *r = &res[3 * (p.off + k)];
                b[k] = r[0]; t[k] = r[1]; v[k] = r[2];
            }
            Done d;
            d.j = p.j; d.ci = p.ci;
            ph_pair_reduce(&cl[p.ci], p.m, b.data(), t.data(), v.data(), &d.r);
            done.push_back(d);
        }
        pend.clear(); cps.clear(); csteps.clear(); cdts.clear();
        return 0;
    };
    // print the finished pairs of targets before `upto`, then checkpoint there
    auto emit = [&](int upto) {
        size_t i = 0;
        for (; i < done.size() && done[i].j < upto; i++) {
            const Done &d = done[i];
            const PhResult &r = d.r;
            fprintf(out, "%s %s %c %d %.3f %.3f %.3f %.4f %.4f %.4f %.2f %.2f %.2f %.4f %.4f\n",
                    tg[d.j].id, cl[d.ci].name, cl[d.ci].kind, r.n, r.d05, r.d50, r.d95, r.t05, r.t50, r.t95,
                    r.v05, r.v50, r.v95, r.f_r50, r.f_rtot);
        }
        done.erase(done.begin(), done.begin() + i);
        if (upto > written) {
            fflush(out);
            FILE *cw = boinc_fopen("phygas_checkpoint", "w");
            if (cw) { fprintf(cw, "%d\n", upto); fclose(cw); }
            boinc_checkpoint_completed();
            boinc_fraction_done((double)upto / nt);
            written = upto;
        }
    };

    std::vector<double> ps((size_t)12 * n_clones);
    for (int j = start; j < nt && !rc; j++) {
        PhTrack tt;
        ph_track_target(&tg[j], tmax, 32, ph_seed("track", tg[j].id), &tt);
        for (int i = 0; i < nc && !rc; i++) {
            if (only[j] ? !in_list(cl[i].name, only[j]) : !ph_prefilter(&tt, &ctr[i], &cl[i])) continue;
            if ((int)csteps.size() + n_clones > cap) {   // device full mid-target
                rc = flush();
                if (!rc) emit(j);
            }
            int steps;
            double dt;
            int m = ph_pair_draw(&tg[j], &cl[i], n_clones, dt_myr, ph_seed(tg[j].id, cl[i].name), ps.data(), &steps, &dt);
            Pend p;
            p.j = j; p.ci = i; p.m = m; p.off = csteps.size();
            pend.push_back(p);
            cps.insert(cps.end(), ps.begin(), ps.begin() + (size_t)12 * m);
            for (int k = 0; k < m; k++) { csteps.push_back(steps); cdts.push_back(dt); }
            n_pairs++; n_cl += m;
        }
        ph_track_free(&tt);
        if (!rc && ((int)csteps.size() >= batch || j + 1 == nt)) {
            rc = flush();
            if (!rc) emit(j + 1);
        }
    }
    if (!rc) fprintf(stderr, "GPU pairs: %d targets, %lld pairs, %lld clones in %.1f s\n", nt - start, n_pairs, n_cl,
                     pg_now() - t_start);
    pg_free(d_prm); pg_free(d_ps); pg_free(d_steps); pg_free(d_dts); pg_free(d_st); pg_free(d_res);
    return rc;
}

static int arg_int(int argc, char **argv, const char *name, int def) {
    for (int i = 1; i + 1 < argc; i++) if (!strcmp(argv[i], name)) return atoi(argv[i + 1]);
    return def;
}
#endif

// target line -> PhTarget, plus the optional "clusters=" list (returned in *only, or empty)
static bool parse_target(char *line, PhTarget *t, char *only, size_t only_size) {
    char f_rv[32], f_erv[32], kind[4];
    int used = 0;
    if (sscanf(line, "%39s %3s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %31s %31s%n",
               t->id, kind, &t->ra, &t->dec, &t->plx, &t->e_plx, &t->dist_pc, &t->e_lndist,
               &t->pmra, &t->e_pmra, &t->pmdec, &t->e_pmdec, f_rv, f_erv, &used) != 14) return false;
    t->kind = kind[0];
    t->rv = strcmp(f_rv, "nan") ? atof(f_rv) : NAN;
    t->e_rv = strcmp(f_erv, "nan") ? atof(f_erv) : NAN;
    t->has_rv = !std::isnan(t->rv) && !std::isnan(t->e_rv);
    only[0] = 0;
    const char *p = strstr(line + used, "clusters=");
    if (p) {
        snprintf(only, only_size, "%s", p + 9);
        only[strcspn(only, " \r\n")] = 0;
    }
    return true;
}

static bool in_list(const char *name, const char *list) {
    size_t n = strlen(name);
    for (const char *p = list; (p = strstr(p, name)) != NULL; p += n)
        if ((p == list || p[-1] == ',') && (p[n] == ',' || p[n] == '\0')) return true;
    return false;
}

int main(int argc, char **argv) {
    BOINC_OPTIONS options;
    boinc_options_defaults(options);
#ifdef PG_GPU
    options.normal_thread_priority = true;  // keep the GPU fed
#endif
    if (boinc_init_options(&options)) return 1;
    fprintf(stderr, "Phygas (birth clusters of runaway stars and neutron stars)\n"
                    "bitboinc -- https://bitboinc.athena.org.tr/\n");
#ifdef PG_GPU
    // device: app_init_data first (BOINC's own assignment), --device as a fallback
    APP_INIT_DATA aid;
    boinc_get_init_data(aid);
    int device = arg_int(argc, argv, "--device", 0);
#ifdef PG_OPENCL
    int platform = arg_int(argc, argv, "--cl-platform", 0);
    if (aid.gpu_opencl_dev_index >= 0) device = aid.gpu_opencl_dev_index;
#else
    int platform = 0;
    if (aid.gpu_device_num >= 0) device = aid.gpu_device_num;
#endif
    if (pg_init(platform, device)) return fail("no usable GPU");
#else
    (void)argc; (void)argv;
#endif

    char in_path[512], cl_path[512], out_path[512];
    if (boinc_resolve_filename("in", in_path, sizeof in_path) ||
        boinc_resolve_filename("clusters", cl_path, sizeof cl_path) ||
        boinc_resolve_filename("out", out_path, sizeof out_path)) return fail("can't resolve file names");

    ph_init();
    PhCluster *cl;
    int nc = ph_read_clusters(cl_path, &cl);
    if (nc <= 0) return fail("can't read the cluster file");

    FILE *in = boinc_fopen(in_path, "r");
    if (!in) return fail("can't open input");
    char header[256], mode[16] = "";
    double p1 = 0, p2 = 0, p3 = 0, p4 = 0;
    if (!fgets(header, sizeof header, in) ||
        sscanf(header, "PHYGAS 1 %15s %lf %lf %lf %lf", mode, &p1, &p2, &p3, &p4) < 3) return fail("bad header");
    bool screen = strcmp(mode, "screen") == 0;
    if (!screen && strcmp(mode, "pair") != 0) return fail("unknown mode");

    // targets
    int cap = 1024, nt = 0;
    PhTarget *tg = (PhTarget *)malloc(sizeof(PhTarget) * cap);
    char **only = (char **)malloc(sizeof(char *) * cap);
    static char line[65536], list[60000];
    while (fgets(line, sizeof line, in)) {
        if (nt == cap) {
            cap *= 2;
            tg = (PhTarget *)realloc(tg, sizeof(PhTarget) * cap);
            only = (char **)realloc(only, sizeof(char *) * cap);
        }
        if (!parse_target(line, &tg[nt], list, sizeof list)) continue;
        only[nt] = list[0] ? strdup(list) : NULL;
        nt++;
    }
    fclose(in);

    double tmax = 0;
    for (int i = 0; i < nc; i++) if (cl[i].tmax_myr > tmax) tmax = cl[i].tmax_myr;

    // per-cluster precomputation (deterministic: seeds come from the names)
    PhLin *clin = NULL;
    PhTrack *ctr = NULL;
    if (screen) {
        ph_screen_r = p3 > 0 ? p3 : 2.0;
        clin = (PhLin *)malloc(sizeof(PhLin) * nc);
        for (int i = 0; i < nc; i++) ph_lin_cluster(&cl[i], (int)p2, ph_seed("lin", cl[i].name), &clin[i]);
    } else {
        ctr = (PhTrack *)malloc(sizeof(PhTrack) * nc);
        for (int i = 0; i < nc; i++) ph_track_cluster(&cl[i], 32, ph_seed("track", cl[i].name), &ctr[i]);
    }

    int start = 0;
    FILE *cp = boinc_fopen("phygas_checkpoint", "r");
    if (cp) {
        if (fscanf(cp, "%d", &start) != 1) start = 0;
        fclose(cp);
        fprintf(stderr, "resuming from checkpoint: target %d of %d\n", start, nt);
    }
    FILE *out = boinc_fopen(out_path, start > 0 ? "a" : "w");
    if (!out) return fail("can't open output");

#ifdef PG_GPU
    if (screen) {
        if (gpu_screen(tg, nt, start, cl, nc, clin, tmax, p1, (int)p2, p4, out)) return fail("GPU screening failed");
    } else {
        if (gpu_pair(tg, only, nt, start, cl, nc, ctr, tmax, (int)p1, p2, out)) return fail("GPU pair run failed");
    }
    start = nt;
#endif
    for (int j = start; j < nt; j++) {
        if (screen) {
            PhLin tl;
            ph_lin_target(&tg[j], tmax, (int)p2, ph_seed("lin", tg[j].id), &tl);
            for (int i = 0; i < nc; i++) {
                int kb;
                double sb, chi2 = ph_screen_chi2(&tl, &clin[i], &cl[i], &kb, &sb);
                if (!(chi2 < p1)) continue;   // dropped anyway: skip the 512-draw probabilities
                PhScreen s;
                ph_screen_finish(&tl, &clin[i], &cl[i], ph_seed(tg[j].id, cl[i].name), chi2, kb, sb, &s);
                if (s.p_r50 < p4) continue;
                fprintf(out, "%s %s %c %.3f %.3f %.2f %.2f %.4f %.4f %.2f\n", tg[j].id, cl[i].name, cl[i].kind,
                        s.chi2, s.t, s.d_nom_pc, s.sig_pc, s.p_r50, s.p_rtot, s.vrel);
            }
            ph_lin_free(&tl);
        } else {
            PhTrack tt;
            ph_track_target(&tg[j], tmax, 32, ph_seed("track", tg[j].id), &tt);
            for (int i = 0; i < nc; i++) {
                if (only[j] ? !in_list(cl[i].name, only[j]) : !ph_prefilter(&tt, &ctr[i], &cl[i])) continue;
                PhResult r;
                ph_pair(&tg[j], &cl[i], (int)p1, p2, ph_seed(tg[j].id, cl[i].name), &r);
                fprintf(out, "%s %s %c %d %.3f %.3f %.3f %.4f %.4f %.4f %.2f %.2f %.2f %.4f %.4f\n",
                        tg[j].id, cl[i].name, cl[i].kind, r.n, r.d05, r.d50, r.d95, r.t05, r.t50, r.t95,
                        r.v05, r.v50, r.v95, r.f_r50, r.f_rtot);
            }
            ph_track_free(&tt);
        }
        if (boinc_time_to_checkpoint() || j + 1 == nt) {
            fflush(out);
            FILE *cw = boinc_fopen("phygas_checkpoint", "w");
            if (cw) { fprintf(cw, "%d\n", j + 1); fclose(cw); }
            boinc_checkpoint_completed();
        }
        boinc_fraction_done((double)(j + 1) / nt);
    }
    fprintf(out, "# done %d\n", nt);
    fclose(out);
    boinc_finish(0);
    return 0;
}
