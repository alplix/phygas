#!/bin/bash
# One-command A/B on a Mac: BASE_REF built the documented way (AGENTS.md "Building": C core -O2 -std=c99, app -O2)
# against this working tree, on the profiles in scripts/bench/profiles. CPU app only: Phygas has no Metal backend, and
# its OpenCL app needs an fp64 GPU, which Apple Silicon has none of (docs/BENCHMARK.md).
# usage: scripts/bench/mac_ab.sh [rounds per profile, default 3] [profiles, default "screen pulsar stage2"]
#   env: BASE_REF     the baseline (default upstream/master if that remote branch exists, else origin/master)
#        PHYGAS_WORK  (default ~/phygas-work): builds in build/ and bins/, results in bench-runs/<date>-<profile>/
#        NEW_OPT, NEW_ALLOC, NEW_PATCHES   how the working tree is built (default: like the baseline, -O2, system malloc,
#                     no patch; NEW_PATCHES is a space-separated list of git-diff patches). With none of them and no
#                     change in src/, the run measures the noise between two identical builds, which is worth doing once
#        FLAGS=1      also builds the working tree with three flag sets and runs them as extra labels: o3 (-O3), lto
#                     (-O2 -flto) and, on Apple Silicon, m1 (-O2 -mcpu=apple-m1 -ffp-contract=off, which keeps the
#                     arithmetic as it is). They are the zero-code, same-results candidates of docs/BENCHMARK.md
#        RUN_PREFIX, GPU_SAMPLE   passed on to run_bench.sh
# Needs the Xcode command line tools and python3, and for NEW_ALLOC=jemalloc Homebrew's jemalloc. Building takes seconds;
# one round of the default profiles takes 2 to 3 minutes per label on an M1 (estimated, not measured). Add the long
# profile with: scripts/bench/mac_ab.sh 3 "screen pulsar stage2 screen-long" (a further 45 minutes).
#
# Labels:
#   base  BASE_REF, the documented build
#   new   the working tree (tracked files, as on disk), built as NEW_OPT, NEW_ALLOC and NEW_PATCHES say
set -e -o pipefail
here="$(cd "$(dirname "$0")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
ROUNDS=${1:-3}
PROFILES=${2:-screen pulsar stage2}
export PHYGAS_WORK=${PHYGAS_WORK:-$HOME/phygas-work}
B=$PHYGAS_WORK/bins
if [ -z "${BASE_REF:-}" ]; then
  if git -C "$repo" rev-parse -q --verify upstream/master > /dev/null; then BASE_REF=upstream/master; else BASE_REF=origin/master; fi
fi
git -C "$repo" rev-parse -q --verify "$BASE_REF^{commit}" > /dev/null || { echo "no git ref $BASE_REF: set BASE_REF"; exit 2; }
echo "baseline $BASE_REF ($(git -C "$repo" rev-parse --short "$BASE_REF")), candidate: working tree on $(git -C "$repo" rev-parse --short HEAD)"
if git -C "$repo" diff --quiet "$BASE_REF" -- src && [ -z "${NEW_OPT:-}${NEW_ALLOC:-}${NEW_PATCHES:-}" ] && [ "${FLAGS:-0}" != 1 ]; then
  echo "note: src/ equals $BASE_REF's and no NEW_* or FLAGS is set: this run measures the noise between identical builds"
fi

REF=$BASE_REF bash "$here/build_variant.sh" base
# NEW_PATCHES is deliberately unquoted: a list of file names
OPT="${NEW_OPT:--O2}" ALLOC="${NEW_ALLOC:-system}" bash "$here/build_variant.sh" new ${NEW_PATCHES:-}
LABELS="base=$B/base/cpu new=$B/new/cpu"
if [ "${FLAGS:-0}" = 1 ]; then
  OPT="-O3" bash "$here/build_variant.sh" o3
  OPT="-O2 -flto" bash "$here/build_variant.sh" lto
  LABELS="$LABELS o3=$B/o3/cpu lto=$B/lto/cpu"
  if [ "$(uname -m)" = arm64 ]; then
    OPT="-O2 -mcpu=apple-m1 -ffp-contract=off" bash "$here/build_variant.sh" m1
    LABELS="$LABELS m1=$B/m1/cpu"
  fi
fi

for prof in $PROFILES; do
  echo "=== profile $prof, $ROUNDS rounds ==="
  bash "$here/run_bench.sh" "$prof" "$ROUNDS" $LABELS
done
