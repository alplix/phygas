#!/bin/bash
# Interleaved speed benchmark of Phygas apps built against scripts/bench/boinc_stub (build_variant.sh).
# usage: run_bench.sh <profile> <rounds> <label>=<binary> [<label>=<binary> ...]
#   profile: a file in scripts/bench/profiles/ without .wu, or a path to a .wu file. Its lines are `key value`: header (the
#        first line of `in`), clusters, targets, kind (S or P), listed, seed (make_bench_data.py), or instead of the
#        last five clusters_file and in_file (real inputs; relative paths from the repository root); # starts a comment
#   Makes the profile's inputs once (make_bench_data.py, or the profile's own clusters_file and in_file), then runs
#   A B C ... A B C ... (rounds) in fresh run directories, so that load drift hits every binary alike. A fresh directory
#   matters here: the app resumes from a phygas_checkpoint it finds.
#   env: BENCH_OUT (default $PHYGAS_WORK/bench-runs, PHYGAS_WORK default ~/phygas-work), RUN_PREFIX (e.g. "taskpolicy -b"
#        or "nice -n 19"), PHYGAS_SRC (where a relative clusters_file or in_file is looked up; default this repository),
#        GPU_SAMPLE=1 (fill the GPU% column: ioreg on macOS, nvidia-smi elsewhere, every 0.5 s; off by default because the
#        sampler is a process of its own, and the CPU app does not use the GPU). A label may name a wrapper script that
#        sets environment variables (PHYGAS_BATCH_STARS, PHYGAS_BATCH_CLONES, PHYGAS_STUB_CHECKPOINT_SECONDS) or passes
#        --device and --cl-platform, and execs the binary.
#   macOS and Linux (GNU time)
# Output: one line per run (label, round, wall s, user s, sys s, max RSS MB, exit code, hash of the result rows, GPU busy %,
# number of result rows), then the summary of summarize_bench.py. The hash covers the rows of `out` (every line but the
# `#` ones) and nothing else: stdout and stderr, which hold the timings, are not hashed. A run without its `# done` line
# gets the hash "incomplete".
set -u
here="$(cd "$(dirname "$0")" && pwd)"
repo=${PHYGAS_SRC:-$(cd "$here/../.." && pwd)}
[ $# -ge 3 ] || { sed -n 2,20p "$0"; exit 2; }
PROFILE=$1; ROUNDS=$2; shift 2
WU=$PROFILE; [ -f "$WU" ] || WU=$here/profiles/$PROFILE.wu
[ -f "$WU" ] || { echo "no profile $PROFILE"; exit 2; }
for lb in "$@"; do
  case $lb in *=*) ;; *) echo "arguments after the rounds are label=binary, not $lb"; exit 2 ;; esac
  [ -x "${lb#*=}" ] || { echo "not an executable: ${lb#*=}"; exit 2; }
done

key() { awk -v k="$1" '$1==k { sub(/^[^ \t]+[ \t]+/, ""); sub(/[ \t]+$/, ""); print; exit }' "$WU"; }
md5_of() { if command -v md5sum > /dev/null; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi; }
hash12() { if command -v shasum > /dev/null; then shasum | cut -c1-12; else sha1sum | cut -c1-12; fi; }
cpu_name() { sysctl -n machdep.cpu.brand_string 2> /dev/null || awk -F': ' '/model name/ { print $2; exit }' /proc/cpuinfo 2> /dev/null || uname -m; }

HEADER=$(key header)
[ -n "$HEADER" ] || { echo "profile $WU needs a 'header PHYGAS 1 ...' line"; exit 2; }
STAMP=$(date +%Y%m%d-%H%M%S)
OUT=${BENCH_OUT:-${PHYGAS_WORK:-$HOME/phygas-work}/bench-runs}/$STAMP-$(basename "$WU" .wu)
mkdir -p "$OUT/data"
CLFILE=$(key clusters_file); INFILE=$(key in_file)
if [ -n "$CLFILE" ] && [ -n "$INFILE" ]; then       # real inputs named by the profile
  case "$CLFILE" in /*) ;; *) CLFILE=$repo/$CLFILE ;; esac
  case "$INFILE" in /*) ;; *) INFILE=$repo/$INFILE ;; esac
  cp "$CLFILE" "$OUT/data/clusters" && cp "$INFILE" "$OUT/data/in" || { echo "can't read $CLFILE or $INFILE"; exit 2; }
else
  NCL=$(key clusters); NTG=$(key targets); KIND=$(key kind); LISTED=$(key listed); SEED=$(key seed)
  [ -n "$NCL" ] && [ -n "$NTG" ] || { echo "profile $WU needs 'clusters <n>' and 'targets <n>' (or clusters_file and in_file)"; exit 2; }
  python3 -B "$here/make_bench_data.py" "$OUT/data" --header "$HEADER" --clusters "$NCL" --targets "$NTG" \
    --kind "${KIND:-S}" --listed "${LISTED:-0}" --seed "${SEED:-1}" || exit 2
fi
[ -s "$OUT/data/clusters" ] && [ -s "$OUT/data/in" ] || { echo "no inputs in $OUT/data"; exit 2; }
NT=$(( $(wc -l < "$OUT/data/in") - 1 ))
RESULTS=$OUT/results.tsv
if [ "$(uname)" = Darwin ]; then
  TIMECMD=(/usr/bin/time -l); RSSDIV=1048576                                   # bytes
else
  TIMECMD=(/usr/bin/time -f "%e real %U user %S sys
%M maximum resident set size"); RSSDIV=1024                                    # GNU time, KiB
fi
{
  echo "# $(date '+%Y-%m-%d %H:%M:%S') profile $WU rounds $ROUNDS cpu $(cpu_name) os $(uname -sr)"
  echo "# header: $HEADER   targets: $NT"
  echo "# data: clusters md5 $(md5_of "$OUT/data/clusters") in md5 $(md5_of "$OUT/data/in")"
  echo "# load: $(uptime | sed 's/.*load/load/')"
  for lb in "$@"; do
    echo "# ${lb%%=*} = ${lb#*=} md5 $(md5_of "${lb#*=}")"
    [ -f "$(dirname "${lb#*=}")/build.info" ] && sed 's/^/#     /' "$(dirname "${lb#*=}")/build.info"
  done
} | tee "$RESULTS"

run_one() { # label binary round
  local label=$1 bin=$2 round=$3 dir=$OUT/$1-$3
  mkdir -p "$dir"
  cp "$OUT/data/clusters" "$OUT/data/in" "$dir/"
  local sampler=""
  if [ "${GPU_SAMPLE:-0}" = 1 ]; then
    if command -v ioreg > /dev/null; then
      ( while :; do ioreg -r -d 1 -c IOAccelerator | grep -o '"Device Utilization %"=[0-9]*' | head -1 | cut -d= -f2; sleep 0.5; done ) > "$dir/gpu_util.txt" 2> /dev/null &
      sampler=$!
    elif command -v nvidia-smi > /dev/null; then
      ( while :; do nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader,nounits 2> /dev/null | head -1; sleep 0.5; done ) > "$dir/gpu_util.txt" 2> /dev/null &
      sampler=$!
    fi
  fi
  ( cd "$dir" && "${TIMECMD[@]}" ${RUN_PREFIX:-} sh -c 'exec "$@" > stdout.txt 2> stderr.txt' sh "$bin" ) 2> "$dir/time.txt"
  echo $? > "$dir/exit_code"
  [ -n "$sampler" ] && { kill $sampler; wait $sampler; } 2> /dev/null
  local real user sys rss gpu=- h=none rows=-
  real=$(awk '/ real /{print $1}' "$dir/time.txt"); user=$(awk '/ real /{print $3}' "$dir/time.txt"); sys=$(awk '/ real /{print $5}' "$dir/time.txt")
  rss=$(awk -v d=$RSSDIV '/maximum resident set size/{printf "%.0f", $1/d}' "$dir/time.txt")
  [ -s "$dir/gpu_util.txt" ] && gpu=$(awk '{s+=$1; n++} END {printf "%.0f", s/n}' "$dir/gpu_util.txt")
  if [ -f "$dir/out" ]; then
    h=$(grep -v '^#' "$dir/out" | hash12)
    rows=$(grep -vc '^#' "$dir/out")
    tail -n 1 "$dir/out" | grep -q '^# done [0-9][0-9]*$' || { h=incomplete; echo "warning: $label round $round: out has no '# done' line" >&2; }
    [ "$rows" = 0 ] && echo "warning: $label round $round: no result rows, so the hash proves nothing (see make_bench_data.py)" >&2
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$label" "$round" "$real" "$user" "$sys" "$rss" "$(cat "$dir/exit_code")" "$h" "$gpu" "$rows" | tee -a "$RESULTS"
  rm -f "$dir/clusters" "$dir/in" "$dir/phygas_checkpoint"
}
for r in $(seq 1 "$ROUNDS"); do
  for lb in "$@"; do run_one "${lb%%=*}" "${lb#*=}" "$r"; done
done
python3 -B "$here/summarize_bench.py" "$RESULTS" | tee "$OUT/summary.txt"
