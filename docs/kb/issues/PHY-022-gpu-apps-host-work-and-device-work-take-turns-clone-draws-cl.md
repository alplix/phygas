---
id: PHY-022
title: GPU apps: host work and device work take turns (clone draws, cluster precompute, result read-back)
type: bug
status: confirmed
severity: low
platforms: [cuda, opencl]
area: performance
action:
next: Measure the host and device shares of a batch on a GPU host (the stderr timing line), then overlap: draw the next flush's clones while the device traces the current one
blocked_by: []
pattern: PAT-033
aliases: []
reported: 2026-11-04
source: read 2026-11-04 (docs/BENCHMARK.md, candidate 2)
where: src/main_phygas_boinc.cpp (gpu_screen, gpu_pair) @ 8003821
fixed_in:
---

## Symptom

Wasted time, no wrong result, not measured: in pair mode the host draws every clone of a flush (about 1 us each, about
0.26 s per flush of 262,144 clones) and only then does the device trace them; in screen mode the cluster precompute,
the 64 draws per star and the read-back and post-processing of each batch run while the device waits. If the two take
similar time, the GPU is idle about half the time.

## Cause

`gpu_screen` and `gpu_pair` (`src/main_phygas_boinc.cpp`) run "prepare, launch, wait, read, process" one after the
other (PAT-033). Read from the code on 2026-11-04 (docs/BENCHMARK.md, "Where the time goes", GPU apps).

## Fix

Proposed, not done: measure the shares first (rule 2 of docs/BENCHMARK.md), then double-buffer: draw the next flush
while the device traces the current one, post-process the previous batch during the device's work. The rows and
their order must stay the same, and the checkpoint after each batch keeps its meaning (PHY-001, PHY-003). Needs a GPU
host with fp64; Apple Silicon cannot run the GPU apps.

## Test

On a GPU host: the same rows before and after, an interleaved A/B on `pulsar` and `screen`.

## Log

- 2026-11-04: recorded with pattern PAT-033.
