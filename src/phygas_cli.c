// Phygas command-line search (testing, validation, and small target lists run on the server).
//   phygas_cli clusters.txt targets.txt [n_clones=1000] [dt_myr=0.05] [only=ID1,ID2,...]
// One line per (target, cluster) pair that passes the prefilter, on stdout:
//   target cluster kind n d05 d50 d95 t05 t50 t95 v05 v50 v95 f_r50 f_rtot age16 age84
// (pc, Myr ago, km/s; f_* = fraction of clones closer than the cluster's r50 / total radius).
// Progress and counts go to stderr.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "phygas.h"

static int wanted(const char *id, const char *only) {
    if (!only) return 1;
    size_t n = strlen(id);
    for (const char *p = only; (p = strstr(p, id)) != NULL; p += n) {
        if ((p == only || p[-1] == ',') && (p[n] == ',' || p[n] == '\0')) return 1;
    }
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s clusters.txt targets.txt [n_clones] [dt_myr] [only=ID,...]\n", argv[0]);
        return 1;
    }
    int n_clones = argc > 3 ? atoi(argv[3]) : 1000;
    double dt = argc > 4 ? atof(argv[4]) : 0.05;
    const char *only = argc > 5 && strncmp(argv[5], "only=", 5) == 0 ? argv[5] + 5 : NULL;
    ph_init();
    PhCluster *cl; PhTarget *tg;
    int nc = ph_read_clusters(argv[1], &cl), nt = ph_read_targets(argv[2], &tg);
    if (nc <= 0 || nt <= 0) { fprintf(stderr, "can't read inputs (%d clusters, %d targets)\n", nc, nt); return 1; }
    double tmax = 0;
    for (int i = 0; i < nc; i++) if (cl[i].tmax_myr > tmax) tmax = cl[i].tmax_myr;
    PhTrack *ctr = malloc(sizeof(PhTrack) * nc);
    for (int i = 0; i < nc; i++) ph_track_cluster(&cl[i], 32, ph_seed("track", cl[i].name), &ctr[i]);
    fprintf(stderr, "%d clusters (traced back up to %.0f Myr), %d targets, %d clones, dt %.3f Myr\n",
            nc, tmax, nt, n_clones, dt);
    long pairs = 0, done = 0;
    clock_t t0 = clock();
    for (int j = 0; j < nt; j++) {
        if (!wanted(tg[j].id, only)) continue;
        PhTrack tt;
        ph_track_target(&tg[j], tmax, 32, ph_seed("track", tg[j].id), &tt);
        int kept = 0;
        for (int i = 0; i < nc; i++) {
            if (!ph_prefilter(&tt, &ctr[i], &cl[i])) continue;
            kept++;
            PhResult r;
            ph_pair(&tg[j], &cl[i], n_clones, dt, ph_seed(tg[j].id, cl[i].name), &r);
            printf("%s %s %c %d %.3f %.3f %.3f %.4f %.4f %.4f %.2f %.2f %.2f %.4f %.4f %.3f %.3f\n",
                   tg[j].id, cl[i].name, cl[i].kind, r.n, r.d05, r.d50, r.d95, r.t05, r.t50, r.t95,
                   r.v05, r.v50, r.v95, r.f_r50, r.f_rtot, cl[i].age16, cl[i].age84);
        }
        fflush(stdout);
        ph_track_free(&tt);
        pairs += kept; done++;
        fprintf(stderr, "%s: %d of %d clusters reachable (%.1f s so far)\n", tg[j].id, kept, nc,
                (double)(clock() - t0) / CLOCKS_PER_SEC);
    }
    fprintf(stderr, "# done %ld targets, %ld pairs, %.1f s CPU\n", done, pairs, (double)(clock() - t0) / CLOCKS_PER_SEC);
    return 0;
}
