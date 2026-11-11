---
id: PHY-021
title: Stage-1 scan computes chi2 at every point, also where a bound shows it cannot pass the cut
type: bug
status: fix-proposed
severity: low
platforms: [all]
area: performance
action: review
next: Review the early reject (ph_screen_chi2_cut); port it to pg_chi2 in phygas_kernel_body.h when an fp64 GPU can measure it
blocked_by: []
pattern:
aliases: []
reported: 2026-11-04
source: M4 benchmark 2026-11-04 (docs/BENCHMARK.md, Results)
where: src/phygas.c:359-385 @ 8003821
fixed_in:
---

## Symptom

Wasted time, no wrong result: stage 1 (screen mode, most of the project's compute) evaluates two 3x3 inverse quadratic
forms and a segment minimum at every point of every (star, cluster) track, although most pairs are far apart compared
with their uncertainty and are dropped by `chi2_max` anyway.

## Cause

`ph_screen_chi2()` (`src/phygas.c:359-385` at 8003821) did not know the cut; the callers (`main_phygas_boinc.cpp:437`,
`phygas_screen.c:55`) drop the pair afterwards when chi2 is not below `chi2_max`. For the symmetric positive-definite M
of a point, m^T M^-1 m >= |m|^2 / trace(M), so a point whose segment stays farther than sqrt(chi2_max * trace(M)) from
the origin cannot bring the pair under the cut.

## Fix

`ph_screen_chi2_cut()` (this branch): skips a point when the squared distance of its step's segment from the origin is
above 2 * chi2_max * trace(M), tested without a division. A skipped point cannot be the minimum of a kept pair, so kept
pairs get the same chi2, kb and sb bit for bit, and dropped pairs stay dropped. `ph_screen_chi2()` is the same function
with an infinite cut. A Gershgorin bound instead of the trace rejected no more points. The GPU kernel (`pg_chi2` in
`phygas_kernel_body.h`) still scans every point and gives the same rows; the same test there is the next step.

## Test

`FLAGS=1 scripts/bench/mac_ab.sh 3 "screen pulsar stage2"`: `screen` 10.96 s -> 9.45 s on an Apple M4 Pro, 16.32 s ->
13.12 s on Linux arm64 GCC 12, the same 6387 rows; `cross_check.sh`: Windows (MinGW-w64) and Linux x86_64 builds write
the same rows. The synthetic stars are a quarter planted near clusters, so the 512-draw finish weighs more than on real
data and the gain on a real chunk should be larger (the scan alone: 1.43x fewer samples).

## Log

- 2026-11-04: created.
- 2026-11-04: fix written and measured (docs/BENCHMARK.md, Results).
- 2026-11-04: phygas_screen CLI (Apple clang 21) on the screen profile's data: 299042 pairs under chi2_max 25 before and after, outputs identical (cmp).
