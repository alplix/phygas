#!/usr/bin/env python3
"""Hunt & Reffert (2024, A&A 686, A42; VizieR J/A+A/686/A42) clusters.dat(.gz) -> Phygas cluster file.

One line per cluster, whitespace separated (names have spaces replaced by '_'):
  idx name type ra dec dist_pc e_dist_pc pmra e_pmra pmdec e_pmdec rv e_rv r50_pc rtot_pc
  age16_myr age50_myr age84_myr tmax_myr
rv / e_rv are 'nan' when the cluster has no Gaia DR3 mean radial velocity (the search then draws
it from the Galactic-rotation expectation, see src/phygas.c).

Kept: open clusters (o) and moving groups (m) that can be younger than MAX_AGE_MYR (16th age
percentile below it). tmax_myr = how far back to trace = 84th age percentile, capped at MAX_AGE_MYR:
a star or neutron star cannot have left a cluster before the cluster formed.
Globular clusters (g) go to a separate file for the millisecond-pulsar search.

Usage: prepare_clusters.py clusters.dat.gz young_clusters.txt globular_clusters.txt
"""
import gzip
import math
import sys

MAX_AGE_MYR = 100.0


def col(line, a, b):
    """Fixed-width field, 1-based inclusive byte range as in the VizieR ReadMe."""
    s = line[a - 1:b].strip()
    return s


def num(line, a, b):
    s = col(line, a, b)
    return float(s) if s else math.nan


def main():
    src, out_young, out_glob = sys.argv[1:4]
    opener = gzip.open if src.endswith(".gz") else open
    n_young = n_glob = 0
    with opener(src, "rt") as f, open(out_young, "w") as fy, open(out_glob, "w") as fg:
        for line in f:
            name = col(line, 1, 20).replace(" ", "_")
            kind = col(line, 281, 281)
            ra, dec = num(line, 320, 331), num(line, 333, 344)
            r50, rtot = num(line, 419, 431), num(line, 461, 473)
            pmra, e_pmra = num(line, 475, 487), num(line, 501, 510)
            pmde, e_pmde = num(line, 512, 523), num(line, 537, 546)
            d16, d50, d84 = num(line, 584, 598), num(line, 600, 614), num(line, 616, 631)
            rv, s_rv, e_rv = num(line, 692, 704), num(line, 706, 718), num(line, 720, 732)
            n_rv = col(line, 734, 737)
            n_rv = int(n_rv) if n_rv else 0
            a16, a50, a84 = num(line, 796, 806), num(line, 808, 818), num(line, 820, 831)
            if any(math.isnan(x) for x in (ra, dec, pmra, pmde, d50)) or d50 <= 0:
                continue
            e_dist = max((d84 - d16) / 2.0 if not math.isnan(d84 - d16) else 0.05 * d50, 0.01 * d50)
            # standard errors can be tiny for rich clusters; floor them (systematics ~0.01 mas/yr)
            e_pmra = max(e_pmra if not math.isnan(e_pmra) else 0.05, 0.01)
            e_pmde = max(e_pmde if not math.isnan(e_pmde) else 0.05, 0.01)
            if math.isnan(rv):
                rv_s, erv_s = "nan", "nan"
            else:
                e = e_rv if not math.isnan(e_rv) else 5.0
                e = max(e, 2.0)
                if n_rv < 3:  # one or two stars: their own scatter says little about the cluster
                    e = max(e, s_rv if not math.isnan(s_rv) else 0.0, 10.0)
                rv_s, erv_s = f"{rv:.3f}", f"{e:.3f}"
            r50 = r50 if not math.isnan(r50) else 2.0
            rtot = rtot if not math.isnan(rtot) else 3.0 * r50
            ages = [10 ** (a - 6) if not math.isnan(a) else math.nan for a in (a16, a50, a84)]
            fields = [name, kind, f"{ra:.6f}", f"{dec:.6f}", f"{d50:.2f}", f"{e_dist:.2f}",
                      f"{pmra:.4f}", f"{e_pmra:.4f}", f"{pmde:.4f}", f"{e_pmde:.4f}", rv_s, erv_s,
                      f"{r50:.3f}", f"{rtot:.3f}"]
            if kind == "g":
                fg.write(f"{n_glob} " + " ".join(fields) + "\n")
                n_glob += 1
                continue
            if kind not in ("o", "m") or math.isnan(ages[0]) or ages[0] >= MAX_AGE_MYR:
                continue
            tmax = min(ages[2] if not math.isnan(ages[2]) else MAX_AGE_MYR, MAX_AGE_MYR)
            fy.write(f"{n_young} " + " ".join(fields) +
                     f" {ages[0]:.3f} {ages[1]:.3f} {ages[2]:.3f} {tmax:.3f}\n")
            n_young += 1
    print(f"{n_young} young clusters/groups, {n_glob} globular clusters")


if __name__ == "__main__":
    main()
