#!/usr/bin/env python3
"""Synthetic Phygas work-unit inputs for the speed benchmark (scripts/bench/run_bench.sh): a cluster file and an `in` file.

The repository holds no catalogue (data/ is ignored: Hunt & Reffert clusters, Gaia stars and ATNF pulsars come from
scripts/prepare_*.py and make_star_catalog.py), so the benchmark generates its inputs. They are a function of the
arguments only: the generator carries its own splitmix64 (the generator of src/montecarlo.h), builds Gaussians from
sums of uniforms, and never uses Python's `random` module or sum(), whose results can change between Python versions.
Item i depends on (seed, i) and not on how many items were asked for, so the 2,000 stars of the quick screen profile are
the first 2,000 of the 20,000 of the long one. The only libm calls are sin, cos, asin, atan2 in the coordinate change;
run_bench.sh records the md5 of both files so a difference between machines shows.

Formats are those of the real writers (clusters: scripts/prepare_clusters.py, 19 fields; stars: make_star_catalog.py;
pulsars: prepare_pulsars.py; header: src/main_phygas_boinc.cpp). The numbers are made up to have the right size and
shape, they are not data:
  clusters  0.15-3.5 kpc from the Sun (median 1 kpc), near the Galactic plane, 85 % open clusters and 15 % moving groups
            (3x larger), 60 % with a radial velocity, age16/50/84 = 0.6/1/1.5 x a median age of 2-100 Myr (log-uniform),
            tmax = age84 capped at 100 Myr, as prepare_clusters.py does (mean tmax about 35 Myr; cluster 0 gets 100, so
            the longest track of a task is 100 Myr as in the real file)
  stars     with a radial velocity (screen mode and stage 2); 25 % are "planted": a straight-line ejection from a random
            cluster 1 Myr to 90 % of its tmax ago (log-uniform) at 5-63 km/s, so that near-misses and result rows exist
            (a pair is kept at chi2 < chi2_max and p_r50 >= p_min, and chance pairs alone may be too rare; the background
            stars are 50-2000 pc from the Sun with 30 km/s random velocities). Planting ignores the Galactic tidal field,
            so some planted pairs end as misses, which is fine
  pulsars   (--kind P) 0.5-5 kpc, a dispersion-measure distance with a 25 % log-normal error, no radial velocity, about
            280 km/s tangential speed, as prepare_pulsars.py writes them
  --listed K  appends `clusters=NAME,...` (K clusters per target: the planted source first) to every target line, which
            is how stage 2 names the clusters to run (no prefilter); without it a pair-mode task runs every cluster the
            prefilter passes (all of them, for a pulsar, whose radial velocity is unknown to +-265 km/s)

usage: make_bench_data.py OUTDIR --header "PHYGAS 1 screen 25 64 2 0.001" --clusters 4087 --targets 2000
                          [--kind S|P] [--listed K] [--seed 1]
Writes OUTDIR/clusters and OUTDIR/in, and a one-line summary on stderr.
"""
import argparse
import math
import os
import sys

M64 = (1 << 64) - 1
PC_PER_MYR = 1.0227121650537077     # 1 km/s in pc/Myr
K_PM = 4.740470446                  # km/s per (mas/yr * kpc)


class Rng:
    """splitmix64 as in src/montecarlo.h."""

    def __init__(self, seed):
        self.s = seed & M64

    def u64(self):
        self.s = (self.s + 0x9E3779B97F4A7C15) & M64
        z = self.s
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & M64
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & M64
        return z ^ (z >> 31)

    def u(self):
        return (self.u64() >> 11) * (1.0 / 9007199254740992.0)

    def uni(self, a, b):
        return a + (b - a) * self.u()

    def gauss(self):
        # Irwin-Hall: only additions, in a fixed order (builtin sum() adds compensated since Python 3.12)
        s = 0.0
        for _ in range(12):
            s += self.u()
        return s - 6.0


def stream(seed, kind, i):
    r = Rng((seed * 0x9E3779B97F4A7C15) ^ (kind * 0xBF58476D1CE4E5B9) ^ ((i + 1) * 0x94D049BB133111EB))
    r.u64()
    return r


# ICRS -> Galactic (the matrix of src/coords.c); the transpose goes back
A_G = ((-0.0548755604162154, -0.8734370902348850, -0.4838350155487132),
       (0.4941094278755837, -0.4448296299600112, 0.7469822444972189),
       (-0.8676661490190047, -0.1980763734312015, 0.4559837761750669))


def gal_to_icrs(v):
    return [A_G[0][i] * v[0] + A_G[1][i] * v[1] + A_G[2][i] * v[2] for i in range(3)]


def norm(v):
    return math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def gal_dir(l, b):
    return [math.cos(b) * math.cos(l), math.cos(b) * math.sin(l), math.sin(b)]


def observables(pos_g, vel_g):
    """Heliocentric Galactic position (pc) and velocity (km/s) -> ra, dec (deg), distance (pc), pmra*, pmdec (mas/yr),
    radial velocity (km/s)."""
    p, v = gal_to_icrs(pos_g), gal_to_icrs(vel_g)
    d = norm(p)
    ra = math.atan2(p[1], p[0])
    if ra < 0:
        ra += 2.0 * math.pi
    dec = math.asin(p[2] / d)
    e_r = [p[0] / d, p[1] / d, p[2] / d]
    e_a = [-math.sin(ra), math.cos(ra), 0.0]
    e_d = [-math.sin(dec) * math.cos(ra), -math.sin(dec) * math.sin(ra), math.cos(dec)]
    k = K_PM * d / 1000.0
    return math.degrees(ra), math.degrees(dec), d, dot(v, e_a) / k, dot(v, e_d) / k, dot(v, e_r)


def make_cluster(seed, i):
    r = stream(seed, 1, i)
    d = 150.0 + 3350.0 * r.u() ** 2
    pos = [d * x for x in gal_dir(r.uni(0.0, 2.0 * math.pi), math.radians(4.0 * r.gauss()))]
    vel = [15.0 * r.gauss() for _ in range(3)]
    a50 = 90.0 if i == 0 else 10.0 ** r.uni(0.3, 2.0)
    a16, a84 = 0.6 * a50, 1.5 * a50
    tmax = min(a84, 100.0)
    kind = "o" if r.u() < 0.85 else "m"
    r50 = 1.0 + 9.0 * r.u() ** 2
    if kind == "m":
        r50 *= 3.0
    rtot = r50 * r.uni(2.0, 4.0)
    ra, dec, dist, pmra, pmdec, rv = observables(pos, vel)
    e_pm = r.uni(0.03, 0.10)
    if r.u() < 0.6:
        rv_s, erv_s = "%.3f" % rv, "%.3f" % r.uni(2.0, 6.0)
    else:
        rv_s, erv_s = "nan", "nan"
    name = "SYN_%05d" % i
    line = "%d %s %s %.6f %.6f %.2f %.2f %.4f %.4f %.4f %.4f %s %s %.3f %.3f %.3f %.3f %.3f %.3f" % (
        i, name, kind, ra, dec, dist, 0.03 * dist, pmra, e_pm, pmdec, e_pm, rv_s, erv_s, r50, rtot,
        a16, a50, a84, tmax)
    return {"name": name, "pos": pos, "vel": vel, "tmax": tmax, "line": line}


def unit_vector(r):
    z = r.uni(-1.0, 1.0)
    phi = r.uni(0.0, 2.0 * math.pi)
    s = math.sqrt(1.0 - z * z)
    return [s * math.cos(phi), s * math.sin(phi), z]


def make_star(seed, j, n_clusters):
    """-> (id, line, index of the source cluster or -1)"""
    r = stream(seed, 2, j)
    pos = vel = None
    src = -1
    if r.u() < 0.25:
        for _ in range(20):
            ci = int(r.u() * n_clusters)
            c = make_cluster(seed, ci)
            t = 10.0 ** r.uni(0.0, math.log10(max(2.0, 0.9 * c["tmax"])))   # log-uniform: short ones keep sigma small
            u = unit_vector(r)
            speed = 10.0 ** r.uni(0.7, 1.8)
            p = [c["pos"][k] + u[k] * speed * PC_PER_MYR * t for k in range(3)]
            if 40.0 <= norm(p) <= 2500.0:
                pos, vel, src = p, [c["vel"][k] + u[k] * speed for k in range(3)], ci
                break
    if pos is None:
        d = 50.0 + 1950.0 * r.u() ** 2
        pos = [d * x for x in gal_dir(r.uni(0.0, 2.0 * math.pi), math.radians(20.0 * r.gauss()))]
        vel = [30.0 * r.gauss() for _ in range(3)]
    ra, dec, dist, pmra, pmdec, rv = observables(pos, vel)
    e_plx, e_pm, e_rv = r.uni(0.02, 0.05), r.uni(0.02, 0.05), r.uni(0.5, 3.0)
    sid = "%d" % (4000000000000000000 + j)
    line = "%s S %.9f %.9f %.6f %.6f 0 0 %.5f %.5f %.5f %.5f %.4f %.4f 12.0 1.0 1.0" % (
        sid, ra, dec, 1000.0 / dist, e_plx, pmra, e_pm, pmdec, e_pm, rv, e_rv)
    return sid, line, src


def make_pulsar(seed, j):
    r = stream(seed, 3, j)
    d = r.uni(500.0, 5000.0)
    pos = [d * x for x in gal_dir(r.uni(0.0, 2.0 * math.pi), math.radians(3.0 * r.gauss()))]
    vel = [200.0 * r.gauss() for _ in range(3)]
    ra, dec, dist, pmra, pmdec, rv = observables(pos, vel)
    e_pm = r.uni(0.3, 1.0)
    tau = 10.0 ** r.uni(0.0, 2.0)
    pid = "JSYN%04d" % j
    line = "%s P %.6f %.6f 0.0000 0.0000 %.1f 0.250 %.4f %.4f %.4f %.4f nan nan %.4g" % (
        pid, ra, dec, dist, pmra, e_pm, pmdec, e_pm, tau)
    return pid, line, -1


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("outdir")
    ap.add_argument("--header", required=True)
    ap.add_argument("--clusters", type=int, required=True)
    ap.add_argument("--targets", type=int, required=True)
    ap.add_argument("--kind", default="S", choices=("S", "P"))
    ap.add_argument("--listed", type=int, default=0)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    if a.clusters < 1 or a.targets < 1 or a.listed < 0 or a.listed > a.clusters:
        sys.exit("make_bench_data.py: need clusters >= 1, targets >= 1, 0 <= listed <= clusters")
    os.makedirs(a.outdir, exist_ok=True)

    tmax_sum = 0.0
    with open(os.path.join(a.outdir, "clusters"), "w") as f:
        for i in range(a.clusters):
            c = make_cluster(a.seed, i)
            tmax_sum += c["tmax"]
            f.write(c["line"] + "\n")

    planted = steps_per_clone = 0.0
    with open(os.path.join(a.outdir, "in"), "w") as f:
        f.write(a.header + "\n")
        for j in range(a.targets):
            tid, line, src = make_star(a.seed, j, a.clusters) if a.kind == "S" else make_pulsar(a.seed, j)
            if src >= 0:
                planted += 1
            if a.listed:
                r = stream(a.seed, 4, j)
                picks = [src] if src >= 0 else []
                while len(picks) < a.listed:
                    ci = int(r.u() * a.clusters)
                    if ci not in picks:
                        picks.append(ci)
                line += " clusters=" + ",".join("SYN_%05d" % ci for ci in picks)
                parts = a.header.split()
                if len(parts) >= 5 and parts[2] == "pair":
                    dt = float(parts[4])
                    for ci in picks:
                        steps_per_clone += math.ceil(make_cluster(a.seed, ci)["tmax"] / dt)
            f.write(line + "\n")

    msg = "%d clusters (mean tmax %.1f Myr), %d %s targets" % (a.clusters, tmax_sum / a.clusters, a.targets,
                                                               "star" if a.kind == "S" else "pulsar")
    if a.kind == "S":
        msg += ", %d planted" % planted
    parts = a.header.split()
    if len(parts) >= 5 and parts[2] == "pair" and a.listed:
        n_clones = int(float(parts[3]))
        msg += ", %d pairs x %d clones = %.3g clone-steps" % (a.targets * a.listed, n_clones, steps_per_clone * n_clones)
    print("make_bench_data.py: " + msg, file=sys.stderr)


if __name__ == "__main__":
    main()
