# Phygas

**Phygas** (φυγάς, "fugitive") traces runaway stars and neutron stars back to the star clusters they
were born in. Every Gaia DR3 star with a full 6-D solution (33.8 million) and every pulsar with a
measured proper motion (~700) is traced back through the Milky Way together with every young star
cluster (4,087 clusters and moving groups, Hunt & Reffert 2024), both drawn from their measurement
uncertainties. The result is an all-sky catalogue of birth-cluster candidates: for each cluster, the
stars and neutron stars that may have left it; for each object, where it may have come from, when,
and how fast it was ejected. It is the sibling of Nostos (interstellar objects) and runs on the same
volunteer-computing project, bitboinc. All raw output is published openly on Zenodo.

## Why

Pulsars against open clusters was done at population level recently (164 pulsars x 2,967 clusters,
arXiv:2511.06850; binary neutron stars, arXiv:2602.14145), but:

1. there is no all-sky catalogue of runaway stars traced to their parent clusters -- existing work
   covers one region or tens of stars (e.g. 39 B stars x 447 clusters, arXiv:2602.16925);
2. runaway-star + neutron-star pairs from the same supernova (a disrupted binary) are searched case
   by case (e.g. zeta Oph + PSR B1706-16, Neuhaeuser et al. 2020);
3. ~700 pulsars have a proper motion, not 145: dispersion-measure distances can enter the Monte Carlo
   with their errors.

## Data

| file | from | rows |
|---|---|---|
| `young_clusters.txt` | Hunt & Reffert 2024 (VizieR J/A+A/686/A42), `scripts/prepare_clusters.py` | 4,087 open clusters + moving groups that can be younger than 100 Myr (878 without a mean RV) |
| `globular_clusters.txt` | same | 132 (for millisecond pulsars, later) |
| `pulsars.txt` | ATNF catalogue, HEASARC copy (Dec 2024), `scripts/prepare_pulsars.py` | 702 with proper motion and distance (356 with a >=3 sigma parallax) |
| `phygas_stars_dr3.txt.gz` | Gaia DR3 6-D (the Nostos download), `scripts/make_star_catalog.py` | ~33.8 M, with errors, Lindegren zero point |

## Method

- Galactic potential, solar motion and coordinate code: the Nostos core (`src/coords.c`,
  `potential.c`, `montecarlo.c`; Allen & Santillan 1991 potential).
- Draws: parallax (or log-normal distance), proper motions, radial velocity. Missing radial velocity:
  clusters -- Galactic-rotation expectation +- 10 km/s; pulsars -- rotation expectation +- 265 km/s
  (Hobbs et al. 2005 birth velocities). Traced back for min(cluster age 84th percentile, 100 Myr).
- **Pulsars: exact Monte Carlo** (`ph_pair`, `phygas_cli`): every cluster that passes a generous
  prefilter, N clones of both, kick-drift-kick leapfrog (dt 0.05 Myr), closest approach refined
  between steps. Output per pair: 5/50/95% of closest distance, time and relative speed (= ejection
  speed), fraction of clones inside r50 / total radius. ~2 min per pulsar at 200 clones.
- **Stars: two stages.**
  1. Linear screening (`ph_lin_*`, `ph_screen`, `phygas_screen`): nominal orbit + state transition
     matrix (6 neighbouring orbits) -> Gaussian position uncertainty at every time; per pair the
     smallest Mahalanobis distance and the probability of being inside r50 / rtot. 0.155 s per star
     against all 4,087 clusters (~60 CPU-days for the whole catalogue).
  2. Exact Monte Carlo (`ph_pair`) for pairs with p_r50 >= 0.01 (~0.6 per star).
- **Chance alignments.** Many pairs are nearby stars crossing large moving groups at 50-80 km/s.
  The analysis repeats the screening with shuffled stellar velocities to measure the false-alarm
  rate per cluster before calling anything a birth site.

## Validation (pulsars, 200 clones)

| case | literature | Phygas |
|---|---|---|
| PSR J1302-6350 (B1259-63) | member of UBC 525 (arXiv:2511.06850) | UBC 525 ranked 1st of 4,087 by p(r50), t ~ 0.1 Myr |
| PSR J0737-3039A/B | Theia 58 more plausible than OC 0450 (arXiv:2602.14145) | Theia 58 in the top 4; OC 0450 weaker |
| PSR J0826+2637 | Stock 7, ~3 Myr (Tetzlaff et al. 2014, pre-Gaia distances) | weak with Gaia distances (Stock 7 at 669 pc): f(rtot) 7.5% |

## First full-size task (20,000 stars, Gaia source_id order, 2026-10-09)

- 51.5 min on one server core (nice 15), 5.8 MB output, 78,122 pairs with p_r50 >= 0.001
  (3.9 per star); validator passes.
- Chance-alignment control (same stars, velocities from other stars):

| cut | real | control | excess |
|---|---|---|---|
| p_r50 >= 0.05, all speeds | 5,318 | 4,835 | 10% |
| p_r50 >= 0.05, relative speed < 20 km/s | 758 | 454 | 67% |

  ~90% of all candidates are chance crossings, but slow pairs (siblings / walkaways) carry a clear
  real excess; among the clusters with the largest excess is the known young beta Tucanae group.
  The catalogue therefore reports, per cluster, real counts against control counts, not raw hits.

## Output and publication

- Pulsar results (one Zenodo record and DOI per pulsar, with a figure each): https://bitboinc.athena.org.tr/phygas_pulsars.php
- Raw: every star-cluster pair with p_r50 >= 0.001 (~3 per star, ~1e8 rows, ~10 GB) and every
  pulsar-cluster pair -> Zenodo, split into sky regions, one concept DOI.
- One Zenodo record (DOI) per cluster (stars and neutron stars that may have left it) and per pulsar.
- Site pages and 3D views like Nostos.

## Next

- BOINC app (`main_phygas_boinc.cpp`): stage-1 work units of ~20,000 stars (~50 min CPU) with the
  cluster file as a sticky input; pulsar work units (one pulsar, 1000 clones, ~10 min).
- Chance-alignment control runs; runaway-star x neutron-star pair search.
