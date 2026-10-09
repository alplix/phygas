#!/usr/bin/env python3
"""ATNF pulsar catalogue (HEASARC copy, atnfpulsar table) -> Phygas target file.

Needs astropy (run with the venv next to the data). Keeps pulsars with a measured proper motion
and some distance. One line per pulsar, whitespace separated:
  id kind ra dec plx e_plx dist_pc e_lndist pmra e_pmra pmdec e_pmdec rv e_rv extra
kind = P (pulsar). Distances: a parallax measured to >= 3 sigma is used as a Gaussian parallax
(plx, e_plx in mas; dist_pc 0); otherwise the catalogue's best distance (mostly from the
dispersion measure, YMW16) with a log-normal error of E_LN_DIST (plx 0). Pulsar radial velocities
are never measured: rv/e_rv are 'nan' and the search draws them (src/phygas.c).
extra = characteristic age P/(2 Pdot) in Myr ('nan' if unknown), for the analysis only.

Usage: prepare_pulsars.py atnfpulsar_heasarc.ecsv pulsars.txt
"""
import math
import sys
import warnings

from astropy.table import Table

E_LN_DIST = 0.25          # ~25% distance error for DM distances (YMW16 is typically 20-40%)
SEC_PER_MYR = 3.15576e13


def val(row, name):
    v = row[name]
    try:
        if getattr(v, "mask", False):
            return math.nan
        return float(v)
    except (TypeError, ValueError):
        return math.nan


def main():
    warnings.simplefilter("ignore")
    t = Table.read(sys.argv[1])
    n = 0
    with open(sys.argv[2], "w") as out:
        for r in t:
            pmra, pmde = val(r, "pm_ra"), val(r, "pm_dec")
            if math.isnan(pmra) or math.isnan(pmde):
                continue
            e_pmra = val(r, "pm_ra_error"); e_pmde = val(r, "pm_dec_error")
            e_pmra = e_pmra if e_pmra > 0 else max(0.1 * abs(pmra), 1.0)
            e_pmde = e_pmde if e_pmde > 0 else max(0.1 * abs(pmde), 1.0)
            plx, e_plx = val(r, "parallax"), val(r, "parallax_error")
            dist = val(r, "distance")
            if plx > 0 and e_plx > 0 and plx / e_plx >= 3:
                p, ep, d = plx, e_plx, 0.0
            elif dist > 0:
                p, ep, d = 0.0, 0.0, dist * 1000.0
            else:
                continue
            per, pdot = val(r, "period"), val(r, "period_dot")
            tau = per / (2 * pdot) / SEC_PER_MYR if per > 0 and pdot > 0 else math.nan
            name = str(r["name"]).replace("PSR ", "").replace(" ", "_")
            out.write(f"{name} P {val(r, 'ra'):.6f} {val(r, 'dec'):.6f} {p:.4f} {ep:.4f} {d:.1f} "
                      f"{E_LN_DIST if d > 0 else 0:.3f} {pmra:.4f} {e_pmra:.4f} {pmde:.4f} {e_pmde:.4f} "
                      f"nan nan {tau:.4g}\n")
            n += 1
    print(f"{n} pulsars with proper motion and distance")


if __name__ == "__main__":
    main()
