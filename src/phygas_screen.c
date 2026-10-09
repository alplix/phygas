// Phygas stage 1 for stars: linear screening of every (star, cluster) pair.
//   phygas_screen clusters.txt stars.txt [chi2_max=25] [n_samp=64] [r_kernel=2 (x r50; 0 = total radius)] [control]
// control: chance-alignment run -- every star keeps its position and parallax but takes the proper
// motions and radial velocity of the star half a file away (so no real birth link survives, only
// chance crossings); comparing a cluster's candidates with and without it gives its false-alarm rate.
// One line per pair with chi2 < chi2_max on stdout:
//   target cluster kind chi2 t_myr d_nom_pc sig_pc p_r50 p_rtot vrel_kms age16 age84
// chi2: smallest Mahalanobis distance over time between the pair's nominal separation and zero,
// with the cluster's total radius added to the covariance (chi2 < ~1: the nominal paths cross inside
// the cluster). p_*: probability of being inside r50 / rtot at that moment. Stderr: progress.
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "phygas.h"

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s clusters.txt stars.txt [chi2_max] [n_samp]\n", argv[0]);
        return 1;
    }
    double chi2_max = argc > 3 ? atof(argv[3]) : 25.0;
    int n_samp = argc > 4 ? atoi(argv[4]) : 64;
    if (argc > 5) ph_screen_r = atof(argv[5]);
    int control = argc > 6 && argv[6][0] == 'c';
    ph_init();
    PhCluster *cl; PhTarget *tg;
    int nc = ph_read_clusters(argv[1], &cl), nt = ph_read_targets(argv[2], &tg);
    if (nc <= 0 || nt <= 0) { fprintf(stderr, "can't read inputs (%d clusters, %d targets)\n", nc, nt); return 1; }
    double tmax = 0;
    for (int i = 0; i < nc; i++) if (cl[i].tmax_myr > tmax) tmax = cl[i].tmax_myr;
    clock_t t0 = clock();
    PhLin *cin = malloc(sizeof(PhLin) * nc);
    for (int i = 0; i < nc; i++) ph_lin_cluster(&cl[i], n_samp, ph_seed("lin", cl[i].name), &cin[i]);
    fprintf(stderr, "%d clusters prepared in %.1f s; %d targets\n", nc, (double)(clock() - t0) / CLOCKS_PER_SEC, nt);
    if (control) {   // velocities of the star nt/2 away; proper motions rescaled by the parallax
        PhTarget *v = malloc(sizeof(PhTarget) * nt);   // ratio so the physical tangential speed is kept
        for (int j = 0; j < nt; j++) {
            v[j] = tg[j];
            const PhTarget *o = &tg[(j + nt / 2) % nt];
            double f = (tg[j].plx > 0 && o->plx > 0) ? tg[j].plx / o->plx : 1.0;
            v[j].pmra = o->pmra * f; v[j].e_pmra = o->e_pmra * f;
            v[j].pmdec = o->pmdec * f; v[j].e_pmdec = o->e_pmdec * f;
            v[j].rv = o->rv; v[j].e_rv = o->e_rv; v[j].has_rv = o->has_rv;
        }
        free(tg); tg = v;
        fprintf(stderr, "control run: velocities taken from the star %d lines away\n", nt / 2);
    }
    long kept = 0;
    t0 = clock();
    for (int j = 0; j < nt; j++) {
        PhLin tl;
        ph_lin_target(&tg[j], tmax, n_samp, ph_seed("lin", tg[j].id), &tl);
        for (int i = 0; i < nc; i++) {
            PhScreen s;
            ph_screen(&tl, &cin[i], &cl[i], ph_seed(tg[j].id, cl[i].name), &s);
            if (!(s.chi2 < chi2_max)) continue;
            kept++;
            printf("%s %s %c %.3f %.3f %.2f %.2f %.4f %.4f %.2f %.3f %.3f\n", tg[j].id, cl[i].name, cl[i].kind,
                   s.chi2, s.t, s.d_nom_pc, s.sig_pc, s.p_r50, s.p_rtot, s.vrel, cl[i].age16, cl[i].age84);
        }
        ph_lin_free(&tl);
        if ((j + 1) % 100 == 0 || j + 1 == nt)
            fprintf(stderr, "%d targets, %ld pairs kept, %.3f s per target\n", j + 1, kept,
                    (double)(clock() - t0) / CLOCKS_PER_SEC / (j + 1));
    }
    return 0;
}
