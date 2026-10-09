#!/usr/bin/env python3
"""Gaia DR3 6D CSV (the one Nostos downloaded) -> Phygas star target catalogue, gzip-compressed.

One line per star with a positive zero-point-corrected parallax, in the format of
scripts/prepare_pulsars.py / src/phygas.h (PhTarget), plus three columns for the analysis/3D view:
  source_id S ra dec plx e_plx 0 0 pmra e_pmra pmdec e_pmdec rv e_rv gmag bp_rp ruwe
Parallax zero point: Lindegren et al. 2021, the same code Nostos uses (nostos/scripts/make_work_units.py).
Unlike Nostos, Phygas draws every star from its own errors (parallax, proper motion, radial velocity).

Usage: python scripts/make_star_catalog.py <gaia_csv> <out.txt.gz> [--nostos-scripts DIR] [--max-rows N]
"""
import argparse
import csv
import gzip
import os
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gaia_csv")
    ap.add_argument("out_gz")
    ap.add_argument("--nostos-scripts", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                            "..", "..", "nostos", "scripts"))
    ap.add_argument("--max-rows", type=int, default=None)
    args = ap.parse_args()
    sys.path.insert(0, os.path.abspath(args.nostos_scripts))
    from make_work_units import zero_points  # noqa: E402

    def g(row, k):
        v = row.get(k, "")
        return v if v not in ("", None) else "nan"

    n_rows = n_written = n_fallback = 0
    pending = []
    tmp = args.out_gz + ".part"
    with gzip.open(tmp, "wt", compresslevel=6, newline="\n") as out:
        def flush():
            nonlocal pending, n_written, n_fallback
            zs, nfb = zero_points(pending, "lindegren")
            n_fallback += nfb
            lines = []
            for row, z in zip(pending, zs):
                try:
                    plx = float(row["parallax"]) - z
                    if plx <= 0:
                        continue
                    lines.append("%d S %.9f %.9f %.6f %.6f 0 0 %.5f %.5f %.5f %.5f %.4f %.4f %s %s %s\n" % (
                        int(row["source_id"]), float(row["ra"]), float(row["dec"]), plx,
                        float(row["parallax_error"]), float(row["pmra"]), float(row["pmra_error"]),
                        float(row["pmdec"]), float(row["pmdec_error"]), float(row["radial_velocity"]),
                        float(row["radial_velocity_error"]), g(row, "phot_g_mean_mag"), g(row, "bp_rp"),
                        g(row, "ruwe")))
                except (ValueError, KeyError, TypeError):
                    continue
            out.write("".join(lines))
            n_written += len(lines)
            pending = []

        with open(args.gaia_csv, newline="") as f:
            for row in csv.DictReader(f):
                n_rows += 1
                if args.max_rows and n_rows > args.max_rows:
                    break
                pending.append(row)
                if len(pending) >= 50000:
                    flush()
                    if n_rows % 1000000 == 0:
                        print(f"{n_rows} rows, {n_written} stars", file=sys.stderr, flush=True)
        flush()
    os.replace(tmp, args.out_gz)
    print(f"done: {n_rows} rows read, {n_written} stars written ({n_fallback} on the global zero-point fallback)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
