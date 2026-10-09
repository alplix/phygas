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

int main(int, char **) {
    BOINC_OPTIONS options;
    boinc_options_defaults(options);
    if (boinc_init_options(&options)) return 1;
    fprintf(stderr, "Phygas (birth clusters of runaway stars and neutron stars)\n"
                    "bitboinc -- https://bitboinc.athena.org.tr/\n");

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

    for (int j = start; j < nt; j++) {
        if (screen) {
            PhLin tl;
            ph_lin_target(&tg[j], tmax, (int)p2, ph_seed("lin", tg[j].id), &tl);
            for (int i = 0; i < nc; i++) {
                PhScreen s;
                ph_screen(&tl, &clin[i], &cl[i], ph_seed(tg[j].id, cl[i].name), &s);
                if (!(s.chi2 < p1) || s.p_r50 < p4) continue;
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
