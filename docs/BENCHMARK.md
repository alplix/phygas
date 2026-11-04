# Benchmarking Phygas (`scripts/bench/`)

How to measure whether a change makes Phygas faster and whether it keeps the results, the scripts that do it, where
the time should go, which changes to try first, and a place for the results. Written 2026-11-04. The method is the one
used for Kampylos (a sibling app of bitboinc); this file repeats it, adapted to Phygas.

**Status.** The suite was written on 2026-11-04 and checked only by reading and by `bash -n` and a Python parse. It has
not been built or run yet: every time in this file is an estimate from the code and from timings recorded elsewhere in
the repository, marked "expected, not measured", and the Results section is empty. The first run should confirm the
estimates (and that the synthetic inputs produce result rows, see "Inputs" below) before anything is concluded.

The scripts build the real BOINC app (`src/main_phygas_boinc.cpp`) against a fake BOINC runtime
(`scripts/bench/boinc_stub/`, the repository has none otherwise) and run it on fixed, generated work units.

## How to run

On a Mac, from the repository root, with the BOINC client paused:

```
scripts/bench/mac_ab.sh              # 3 rounds of the profiles screen, pulsar and stage2
scripts/bench/mac_ab.sh 5 screen     # 5 rounds of one profile
NEW_OPT="-O3" scripts/bench/mac_ab.sh          # the working tree built with other flags
FLAGS=1 scripts/bench/mac_ab.sh                # also -O3, -flto and -mcpu=apple-m1 -ffp-contract=off, as extra labels
```

It needs the Xcode command line tools and `python3` (jemalloc only for `NEW_ALLOC=jemalloc`). It builds, in
`~/phygas-work`, the baseline `upstream/master` (or `origin/master`; `BASE_REF=` picks another ref) the documented way
(AGENTS.md "Building": C core `-O2 -std=c99`, app `-O2`) as label `base`, and this working tree as label `new`, then
runs each profile with `base` and `new` interleaved. With no change in `src/` and no `NEW_*` setting, `new` is
identical to `base`: run that once first, because its spread is the noise a real change has to beat. Each profile ends
with a table (`~/phygas-work/bench-runs/<date>-<profile>/summary.txt`, the runs in `results.tsv`) like this one:

```
label        runs  median s    min s    max s    cpu s  RSS MB  GPU%    rows   vs 1st  verdict / output
base            3         .        .        .        .       .     -       .        .  baseline; same output
new             3         .        .        .        .       .     -       .        .  ACCEPT (...) / no clear gain (...); same output
```

`verdict` applies the rule of this file: a candidate is faster when its median is at least 5 % below the baseline's and
its slowest run beats the baseline's fastest, or when its median is at least 10 % below. `same output` means every run
of that label wrote the same result rows as the baseline (below), `same output as X` the same as an earlier label (a
GPU build after CPU builds: they differ at rounding level, PHY-005), otherwise the column shows the change of the row
count (`OTHER output (rows +3)`). `rows` is the number of result rows, the key number of a run: it must not be 0 (then
the hash proves nothing and the script says so), and a change that moves it changed the results.

By hand, for other builds and other machines (Linux builds the CPU app, and the CUDA app when `nvcc` is found):

```
REF=upstream/master scripts/bench/build_variant.sh base            # bins/base/cpu
OPT="-O2 -flto" scripts/bench/build_variant.sh lto                 # the working tree; a patch: build_variant.sh lto my.patch
scripts/bench/run_bench.sh screen 3 base=$HOME/phygas-work/bins/base/cpu lto=$HOME/phygas-work/bins/lto/cpu
```

`build_variant.sh` copies the tracked files (`git ls-files`; modified ones as they are on disk, a new file must be
`git add -N`ed) or `git archive`s `REF`, so a stale object file or an untracked file cannot leak in. The stub and the
data generator always come from this checkout, so a ref from before `scripts/bench/` (such as `upstream/master` today)
can be built. Knobs: `OPT`, `CC`, `CXX`, `CXXSTD`, `ALLOC` (system, jemalloc, mimalloc), `BACKENDS`, `REF`, and for
the GPU apps `VENDOR`, `NVCC`, `CUDA_ARCH`, `NVFLAGS`; each build writes `bins/<name>/build.info` (compiler, flags, ref),
which `run_bench.sh` copies into `results.tsv`. `run_bench.sh` knobs: `RUN_PREFIX` (`taskpolicy -b`, `nice -n 19`),
`GPU_SAMPLE=1`, and any `PHYGAS_STUB_CHECKPOINT_SECONDS`, `PHYGAS_BATCH_STARS`, `PHYGAS_BATCH_CLONES` in the
environment; a label may name a wrapper script that sets them or passes `--device` and `--cl-platform`.

The GPU apps are written into the scripts (`BACKENDS="cuda"`, `BACKENDS="cpu opencl"`) but have **not been tried**: CUDA
needs the toolkit (the line is `nvcc -O3 --fmad=false` for the device code, `-arch=native` unless `CUDA_ARCH` says
otherwise, not the release build's list of architectures), OpenCL needs an fp64 GPU, so Apple Silicon cannot run it
(the OpenCL build itself compiles on a Mac). Without `nvcc` the CUDA app is a gap of this suite.

## Method

### Two kinds of speed change

* **Same results**: build flags, loop structure that keeps every operation and its order, early exits that rest on an
  exact bound, scheduling of host and GPU work, threads that keep the row order. The result rows (`out` without its `#`
  lines) must stay byte-identical to the old build's on the same toolchain. Such a change needs a measured gain, the
  output check and the correctness checks below, nothing more.
* **Result-changing**: anything that changes a number the app prints: another formula for a force (`pow` rewritten),
  another step size (`dt_myr`, `PH_LIN_SUB`, `PH_PRE_STEP_MYR`), fewer clones, a tighter prefilter, another RNG, another
  order or number of draws, FMA contraction, `-ffast-math` (not allowed, AGENTS.md "Floating point"). Measure it the
  same way, and treat it as a science change that the maintainers decide on: the validator compares replicas (what it
  compares is not known here, PHY-005), published results come from the current apps (pulsar results on Zenodo, one
  record per pulsar), and the orbit core is a copy of nostos's (PHY-012). Until then keep it behind a compile-time
  switch that is off by default.

### Rules for every measurement

1. **Build both sides the same way, from clean copies.** The baseline is a git ref (`REF=upstream/master`) built with the
   documented build line; the candidate is the working tree, or the baseline plus a patch. The stub replaces BOINC's
   libraries on both sides, so the release builds (static, cross-compiled, real BOINC) are not what is measured.
2. **Find out where the time goes first.** Sample a run (`sample <pid> 10` on macOS, `perf record -g` on Linux) and, for
   a GPU app, split GPU time from host time: the GPU app prints `GPU screening: N stars in X s` or `GPU pairs: ... in X s`
   to stderr (`<label>-<round>/stderr.txt`). A change saves at most the share of time it touches: a change to 10 % of
   the time cannot gain more than 1.11x.
3. **Use a profile that exercises the change.** A change to the scan over track points does nothing on `pulsar`; a change
   to the force law does little on `screen` (the integration is a small part of it) and everything on `pulsar` and
   `stage2`. Fixed costs weigh more in a quick profile than in a task (the cluster precompute is about 1 s in a 40 s
   run and 0.3 % of a 20,000-star task). Confirm the final number on `screen-long` (stage 1) and on a pulsar task of 1000
   clones (copy `pulsar.wu` and edit its header).
4. **Interleave and repeat.** `run_bench.sh` runs A B C A B C ... in fresh directories (a fresh directory is also what
   keeps the app from resuming from a `phygas_checkpoint`); use at least 3 rounds. The load from other programs changes
   over minutes, and interleaving makes it hit every build alike. Compare medians.
5. **Accept by the rule, not by eye.** A candidate is faster when its median is at least 5 % below the baseline's and its
   slowest run beats the baseline's fastest, or when its median is at least 10 % below. Anything else is "no clear gain"
   and is reported as that.
6. **Check the output of every run.** Each run's result rows are hashed (every line of `out` but the `#` ones; stdout and
   stderr, which hold timings, are not); a run without its `# done` line is "incomplete". Compare like with like: other
   compilers, other platforms and the GPU apps differ from the CPU app at rounding level (PHY-005; clang that contracts
   `a*b+c` is the case to watch), so a Linux GCC build is compared with a Linux GCC build, and a CUDA build with
   the CUDA build before it. A printed row rarely changes at rounding level, but a row near `chi2_max`, `p_min` or the
   prefilter limit can appear or vanish.
7. **Measure on a quiet machine.** Pause the BOINC client (CPU and GPU tasks) and other heavy programs, keep a laptop on
   power, and look at `uptime` (the scripts record the load). The CPU app is one thread, so one busy core is enough to
   spoil a run. On a loaded machine an interleaved comparison is still fair, but its absolute times and spreads are not
   worth quoting.
8. **Write it down.** Date, machine, OS, compiler, load, profile, rounds, medians and spread, verdict and output go into
   the Results below. A number in a commit message, a README or a code comment comes from there, and a platform that was
   not measured is called unmeasured.

### Steps for one change

1. Sample the current build on the profile closest to production and note where the time goes.
2. Make the change in the working tree, or as a patch (`build_variant.sh <name> <patch>`); build the baseline with `REF=`.
3. Run one round (`run_bench.sh screen 1 ...`) to see that the output stays the same and that `rows` is not 0.
4. Run 3 or more rounds on a quiet machine, on each profile the change touches, and `screen-long` or a 1000-clone pulsar
   task for the final number.
5. Run the correctness checks AGENTS.md routes the touched files to. The generated inputs work with the CLIs too
   (in `bench-runs/<date>-<profile>/data/`: `tail -n +2 in > stars.txt`, then `phygas_screen clusters stars.txt` and
   `phygas_cli`, `cmp` before and after). A
   change to the core or the kernel also needs the same change in `phygas_kernel_body.h` (then
   `python3 scripts/embed_kernel.py`), in nostos's copy of the core (PHY-012), and a GPU run compared with the CPU app. A
   change to checkpointing or to how targets are batched needs the quit-and-kill run (PHY-001); this suite does not do it.
6. Measure other platforms where possible (Linux natively or in a container, the CUDA and OpenCL apps on a GPU host), and
   name the ones that were not measured.
7. Record the results here, and update AGENTS.md "Building" when the way to build changes.

### Traps (read from the code and the knowledge base; none seen in a measurement yet)

* **Seeds and the generated data.** A pair's draws depend on `ph_seed(target id, cluster name)`, so the generated ids and
  names are part of the profile: change the generator and every hash and time changes. Compare only runs of one
  `run_bench.sh` call (it generates the data once and records their md5), and do not edit `make_bench_data.py` between the
  two sides of a measurement. The data are a function of the profile and the Python's libm (`sin`, `cos`); another machine
  may generate slightly different digits.
* **The draws are the results.** Splitmix64 and Box-Muller (`montecarlo.c`, ported in `phygas_kernel_body.h`) and the order
  of the draws in `ph_sample_*`, `ph_pair_draw` and `ph_screen_finish`: reusing the second Box-Muller value, drawing in
  another order, or another seed rule changes every row, even though it saves almost nothing (draws are about 1 % of pair mode).
* **The synthetic inputs are guesses at the real ones.** A profile's pair counts (how many clusters the prefilter passes,
  how many pairs pass `chi2_max` and `p_min`) were not measured. If the `rows` column is 0 or far from "a few per star",
  the profile does not exercise `ph_screen_finish` or the output path the way production does; change the generator or the
  profile's header, or give the profile real files (`clusters_file`, `in_file`, see the usage text of `run_bench.sh`).
* **Checkpoints.** The CPU app checkpoints when `boinc_time_to_checkpoint()` says so: the stub says so every 60 s
  (`PHYGAS_STUB_CHECKPOINT_SECONDS`, 0 = after every target), so a 40 s quick profile writes none until the end, and a
  checkpoint-heavy comparison needs the knob in every run. The GPU paths checkpoint after every batch whatever the timer
  says (PHY-003), so their number of checkpoints follows `PHYGAS_BATCH_STARS` and `PHYGAS_BATCH_CLONES`: set the same
  values in all runs of a comparison. The stub has no timer thread, no heartbeat and no 2 s sleep in `boinc_finish`; their
  cost is not in the numbers. A pulsar task is one target and has no checkpoint (PHY-011).
* **Resume.** A change that changes the batching of targets changes what a quit throws away; a resume today can repeat or
  tear rows (PHY-001), so "same output as an uninterrupted run" is not a property to preserve blindly. This suite never
  quits a run.
* **GPU sharing.** Another GPU task, or the display, on the same GPU slows the device part and the sizing of launches:
  `pg_next_slice` aims at 0.25 s per launch from the last launch's time, so the slices (not the results: each launch
  continues from saved state) depend on the load. Consumer NVIDIA cards run fp64 at a small fraction of their fp32 rate,
  which is why the device part is fp64-bound and the host part can show. Keep the BOINC client's GPU work paused.
* **Thread counts.** The CPU app is one thread, and `APP_INIT_DATA.ncpus` is not read; the GPU app's host side is one
  thread. A threaded variant is not only a speed change: the client schedules an app by its plan class's `avg_ncpus`, so
  more threads would oversubscribe volunteers' machines unless the server changes too. Results do not depend on thread
  count if the row order is kept (the seeds are per pair), so the output check works; the decision is the maintainers'.
* **Efficiency cores.** On Apple Silicon a process with background QoS (`taskpolicy -b`) runs on the efficiency cores, a
  niced one need not, and the BOINC client's way of lowering priority decides which a volunteer's Mac does. Compare
  `RUN_PREFIX="taskpolicy -b"` and `RUN_PREFIX="nice -n 19"` before trusting a gain that matters on volunteers' machines.
* **Published results and the shared core.** `coords`, `potential`, `montecarlo`, `units`, `vec3` are byte-identical
  copies of nostos's (`docs/kb/kb.conf`, PHY-012): a speed change there belongs in both repositories and changes nostos's
  results too. The authors measured the CPU and GPU apps byte-identical (commits 8bb46cb, b63a70a); keep that claim true
  or say that it no longer is.
* **Floating-point flags.** The core is built with `-std=c99` (GCC then does not fuse `a*b+c`), the CUDA code with
  `--fmad=false`, the OpenCL code with `FP_CONTRACT OFF`. Apple clang 14 does not contract by default, newer clang does:
  a flag set or a compiler update can change rounding without any source change (the 2026-10-27 measurement in PHY-005
  found 60,000 screen rows and 1,631 pair rows printed identically with contraction on and off, which is not a guarantee).
* **Dynamic against static, one libm.** The release Linux builds are fully static; the benchmark builds are dynamic and use
  the machine's libm. Ratios between builds on one machine carry over, absolute times and the `pow` share may not.

## The pieces

| file | what it does |
|---|---|
| `build_variant.sh <name> [patch ...]` | copies the tracked files (or `REF=<git ref>`), applies patches, builds `bins/<name>/cpu` (macOS and Linux), `cuda` (when `nvcc` is found, or named) and `opencl` (named), and writes `build.info` |
| `run_bench.sh <profile> <rounds> label=binary ...` | makes the inputs, then interleaved runs in fresh directories; wall, user and system time, peak memory, a hash of the result rows, their count. `RUN_PREFIX` runs every binary through a command |
| `summarize_bench.py <results.tsv>` | the table above |
| `profiles/*.wu` | header line, sizes and seed of one work unit, with the reason for it |
| `make_bench_data.py` | writes the `clusters` and `in` files of a profile, deterministically |
| `boinc_stub/` | the BOINC calls the app makes (`boinc_api.h`, `filesys.h`, `boinc_stub.cpp`), for plain processes |
| `mac_ab.sh` | everything above in one command |

| profile | work unit | why | expected on an M1 |
|---|---|---|---|
| `screen` | stage 1: 2,000 stars against 4,087 clusters, `PHYGAS 1 screen 25 64 2 0.001` | the chi2 scan over every (star, cluster, track point), which should be most of stage 1 on the CPU | 30-60 s |
| `pulsar` | one pulsar against 4,087 clusters, `PHYGAS 1 pair 100 0.05`, no `clusters=` list | the pair loop (`ph_pair_clone`) under the prefilter, which passes nearly every cluster for a pulsar | 30-80 s |
| `stage2` | 30 stars with 10 named clusters each, `PHYGAS 1 pair 1000 0.05` | the same loop with many small pairs, known radial velocities and the `clusters=` path; exact size | 20-60 s |
| `screen-long` | stage 1: 20,000 stars (the first 2,000 are `screen`'s) against 4,087 clusters | one production-size stage-1 task, and the GPU apps, whose 12 s per task makes the quick profile too short | about 7 min (CPU) |

Times are expected, not measured: they come from 0.02 s per star against 4,087 clusters (`docs/kb/ledger.md`, PAT-007: "about
0.02 s screen" on an M1; the README's 0.155 s predates the 8x speed-up of commit 8bb46cb) and from the pair cost in
"Where the time is expected to go".

### Inputs

The repository holds no catalogue, so `make_bench_data.py` writes the inputs of a profile from its sizes and a seed
(own splitmix64, no `random` module; item i depends only on the seed and i, so a smaller profile is a prefix of a
larger one). The formats are the real ones (`prepare_clusters.py`, `make_star_catalog.py`, `prepare_pulsars.py`, the
header of `main_phygas_boinc.cpp`); the numbers are made up with the right shape: clusters 0.15-3.5 kpc from the Sun near
the plane, ages 2-100 Myr with `tmax` as `prepare_clusters.py` computes it (mean about 35 Myr, longest 100), 60 % with a
radial velocity; stars with 6-D data, a quarter of them planted as ejections from a cluster so that near-misses exist;
a pulsar with a dispersion-measure distance and no radial velocity. The headers are not recorded anywhere in the
repository: `chi2_max` 25, `n_samp` 64 and `r_kernel` 2 are the defaults of `phygas_screen`, `p_min` 0.001 is the cut of the
first full-size task in README.md, `dt_myr` 0.05 is the CLI's default, 1000 clones is the pulsar task's. A profile can
name real files instead (`clusters_file`, `in_file`; relative paths from the repository root).

## Where the time is expected to go (expected, not measured)

Read from the code at 0a3eda2; to be confirmed with `sample` or `perf` before a change is built on it (rule 2).

* **CPU app, screen mode** (about 0.02 s per star on an M1, documented):
  * `ph_screen_chi2` for every cluster, about 90 %. The scan visits each point of the cluster's own track (its `tmax` over
    0.5 Myr: about 71 points on average in the profiles, so about 290,000 points per star) and does, per point, a 3x3
    covariance sum, `quad_inv` (about 40 flops and a division), `seg_min` (a division and a square root) and a second
    `quad_inv`. 0.02 s over 290,000 points is about 70 ns per point.
  * `lin_integrate` of the star: 201 points x 10 sub-steps x 7 bodies = 14,070 leapfrog steps, about 1 ms, 5 %.
  * `ph_screen_finish` (512 x 3 Gaussians, each a `log`, a `sqrt` and a `cos`, about 50 us) for every pair that passes
    `chi2_max`; its share depends on the pairs per star, which is unknown for the synthetic data (up to about 5 %).
  * Once per task: `ph_lin_cluster` for 4,087 clusters, about 1.3 s (ledger: 0.1 s per 300 clusters on an M1), 0.3 % of a
    20,000-star task and 3 % of `screen`.
  * Reading, parsing and writing rows: below 1 %.
* **CPU app, pair mode**: more than 90 % in `ph_pair_clone`. Per step two `body_step` calls and one `seg_min`; each
  `body_step` is one `milkyway_accel` through the function pointer in `Potential`, with three `pow` calls (`pow(r^2+b^2, 1.5)`
  and two in the halo), three `sqrt` and several divisions. A guess of about 145 ns per clone-step on an M1, of which the
  `pow` calls may be 40 %. `pulsar`: 4,087 pairs x 100 clones x about 700 steps = 2.9e8 clone-steps, about 40 s (README:
  about 2 minutes per pulsar at 200 clones, machine not named). `ph_pair_draw` is about 12 Gaussians and two coordinate
  transforms per clone, about 1 us against about 100 us of integration: 1 %. `ph_pair_reduce` (three `qsort`s of
  `n_clones` values) and the prefilter: below 1 %. `ph_track_cluster` for 4,087 clusters: about 0.2 s once.
* **GPU app, screen mode** (documented: 20,000 stars in 12 s on an RTX 5070 Ti, commit 8bb46cb; the repository does not
  split it). Device: `pg_chi2` (one thread per star and cluster, the same scan), `pg_lin` (one thread per star: only
  `PHYGAS_BATCH_STARS` = 1024 threads, each 14,070 serial steps, so poor occupancy), `pg_finish`. Host, serial and not
  overlapped with the device: the cluster precompute (about 1.3 s on an M1 core, probably a visible part of 12 s), the
  64 draws per star (`ph_lin_target_init`), reading back `batch x 4,087` results of three arrays (about 85 MB per
  batch), the loop over them with `ph_seed` for every pair that passes `chi2_max`, and the row output. The device part is
  limited by fp64 throughput.
* **GPU app, pair mode**: the host draws every clone (`ph_pair_draw`, about 1 us each, so about 0.26 s per flush of
  262,144 clones) and the device then traces them (`pg_pair_k`, one thread per clone, launches of about 0.25 s); nothing
  overlaps. If the device takes about as long as the host draws, the GPU is idle half the time. Unknown until measured.

## Candidate speed-ups, in the order to measure them

Each is a hypothesis to test with the method above, not a result. "Same" means the result rows should stay byte-identical
(the hash decides), "changes results" means they will not.

1. **Division-free early reject in the stage-1 scan: same results if the bound is exact.** The CPU loop drops a pair whose
   chi2 is not below `chi2_max` and prints chi2, `kb`, `sb` only for kept pairs, so a point need not be evaluated when it
   certainly cannot lower the minimum below `chi2_max`. For symmetric positive-definite M, `m^T M^-1 m >= |m|^2 / trace(M)`,
   and `|m|^2 > chi2_max * trace(M)` needs no division (with a small safety factor against rounding; the segment test
   needs the same bound at the closest point). Most (star, cluster) pairs are far apart compared with their position
   uncertainty, so most points should be rejected, in the code that is about 90 % of stage 1, which is most of the
   project's compute. It needs `chi2_max` passed into `ph_screen_chi2` (the CLI and the app both call it). The GPU kernel is
   a literal port and need not change for the output to stay the same. Measure on `screen`, then `screen-long`.
2. **Overlap host and device work in the GPU apps: same results.** Pair mode: draw the next batch's clones while the device
   traces the current one. Screen mode: prepare the next batch's stars and post-process the last batch's results during the
   device's work. Both are possible gains of the order of the host share, which rule 2 must measure first (the stderr line
   splits it). The checkpoint after each batch has to keep its meaning (PHY-001, PHY-003).
3. **Flag sets: same results if the hash agrees.** `-O3`, `-flto`, `-mcpu=apple-m1` or `-march=native` with
   `-ffp-contract=off` (without it newer clang may contract and change rounding). No code change, so the cheapest test
   (`FLAGS=1 scripts/bench/mac_ab.sh`). For Kampylos such flags gave 0-8 %, within noise; the force function is called
   through a function pointer built at run time (`ph_init`), which LTO may not remove.
4. **Fewer `pow` calls per force evaluation: changes results at rounding level.** `pow(r^2+b^2, 1.5)` as `d * sqrt(d)` and
   the halo's `pow(x, gamma)` as `x * pow(x, gamma-1)` leave one `pow` of three. It is the largest lever on both
   modes (every leapfrog step, in `lin_integrate`, `ph_pair_clone` and the device code), maybe 20-35 % of a step. It
   changes the last bits of every orbit; PHY-005 records that the `x * sqrt(x)` form "would change science results and is
   not proposed without the maintainers". It needs the same change in nostos's core, in `phygas_kernel_body.h` and a
   count of changed rows on a real task.
5. **Contiguous (structure-of-arrays) cluster tracks and SIMD in the scan: same results.** The per-point arithmetic stays
   the same, so values are identical if the minimum is still taken in order. The compiler cannot do this for the current
   array-of-structures layout; the production x86-64 builds have only SSE2 (two doubles). Worth it only if candidate 1 leaves
   the scan large.
6. **Threads on the CPU: same results, a resource decision.** Pairs and stars are independent, and the seeds are per pair,
   so rows do not change if their order is kept. It is a plan-class and server question first (trap "Thread counts").
7. **Step sizes, clone counts, prefilter margin: change results.** `dt_myr` 0.05 to 0.1 halves pair mode, `PH_LIN_SUB` 10 to
   5 halves `lin_integrate`, a tighter prefilter removes pairs. These need a convergence study and the maintainers; measure
   the effect on rows before the effect on time.
8. **Call the force directly in `body_step` (no function pointer, one translation unit): same results**, probably under
   10 %; it can enable inlining that candidate 3's LTO does not.

Not worth measuring first: the allocator (`ALLOC=`: Phygas allocates per target and per pair, not per step, so expected
0 %; it is in the scripts to confirm that), the Gaussian generator (about 1 % of pair mode, and changing it changes
every result), parsing and output.

## What the suite does not cover

* Windows, Android, 32-bit and the other Linux architectures, and the release builds (static, cross-compiled, real BOINC
  libraries): the suite builds natively against the stub.
* The CUDA and OpenCL apps are written into `build_variant.sh` and have not been tried. Without `nvcc` there is no CUDA
  app; Apple Silicon has no fp64 OpenCL device, so the GPU apps cannot be measured on a Mac. A GPU host is needed, with
  `screen-long` and a longer pulsar profile.
* Real data: the Hunt & Reffert clusters, the Gaia stars and the ATNF pulsars are not in the repository. The synthetic pair
  counts, prefilter pass rate and age distribution are estimates, and the headers of real work units come from the server.
* Quit, kill and resume (no test exists, PHY-001), the real client's timer, heartbeat and 2 s `boinc_finish` sleep, and what
  the validator accepts.
* Other work on the volunteer's machine, thermal limits of a phone or laptop, and runs longer than a task.
* The two CLIs (`phygas_cli`, `phygas_screen`): they share the core but have their own loops. The chance-alignment
  `control` mode, and results across platforms (PHY-005).

## Results

Nothing has been measured yet. Template for each entry (copy it under this heading, newest first):

```
### <change>, <date>
Machine <model, cores, RAM>, OS <version>, compiler <version> (bins/<name>/build.info), flags <...>, load <uptime>,
BOINC client <paused / not running>, RUN_PREFIX <none / ...>, profile <name>, rounds <n>, data md5 <from results.tsv>.

label        runs  median s    min s    max s    cpu s  RSS MB  GPU%    rows   vs 1st  verdict / output
<paste the table of summary.txt>

Output: <same as baseline / differs: how many rows, which, why>.
Where the time went: <sample or perf summary, GPU time against host time for a GPU app>.
Checks run: <CLI cmp, kernel-body and nostos copies, GPU against CPU, kill and resume ...>.
Not measured: <platforms, backends, profiles>.
```
