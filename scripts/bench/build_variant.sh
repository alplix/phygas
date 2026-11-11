#!/bin/bash
# Builds one variant of the Phygas BOINC app against the fake BOINC runtime (scripts/bench/boinc_stub), for run_bench.sh.
# usage: build_variant.sh <name> [patch ...]
#   Copies the repository's tracked files (modified ones as they are on disk; a new file needs `git add -N`) to
#   $PHYGAS_WORK/build/<name>, or the files of git ref $REF, applies the patches (git diff format, -p1) and builds
#   $PHYGAS_WORK/bins/<name>/<backend> with the lines of AGENTS.md "Building": C core `-O2 -std=c99`, app `-O2`.
#   The stub and this script come from this checkout also when REF is built (a ref from before scripts/bench has none).
#   env: BACKENDS  default "cpu" on macOS; on Linux "cpu", plus "cuda" when nvcc is found. "opencl" only when named: it
#                  builds on macOS but cannot run on Apple Silicon (no fp64 device, src/phygas_opencl.c refuses it)
#        OPT       flags for every C and C++ compile and the link (default -O2). -std=c99 is always added to the C core;
#                  it and the FP flags of AGENTS.md "Floating point" must stay (no -ffast-math)
#        CC, CXX   compilers (default clang and clang++ on macOS, cc and c++ elsewhere)
#        CXXSTD    standard for the app and the stub (default -std=c++17, which the GPU lines use too; the CPU line of
#                  AGENTS.md has none, CXXSTD= gives that). GPU backends always use -std=c++17
#        ALLOC     system (default), jemalloc (JEMALLOC = static library or linker flags; default Homebrew's
#                  libjemalloc.a on macOS, -ljemalloc elsewhere) or mimalloc (MIMALLOC, default -lmimalloc). Phygas
#                  allocates per target and per pair, not per step: expected to change nothing, measure before keeping
#        REF       build this git ref (e.g. upstream/master) instead of the working tree
#        VENDOR    opencl: the -DEXPECTED_GPU_VENDOR of the release builds ("nvidia" or "amd"); empty = any platform
#        NVCC, CUDA_HOME, CUDA_ARCH (default native, needs nvcc 11.5), NVFLAGS   cuda; device flags stay -O3 --fmad=false
#        PHYGAS_WORK (default ~/phygas-work), PHYGAS_SRC (the checkout to copy; default this repository)
set -e -o pipefail
here="$(cd "$(dirname "$0")" && pwd)"
repo=${PHYGAS_SRC:-$(cd "$here/../.." && pwd)}
[ $# -ge 1 ] || { sed -n 2,21p "$0"; exit 2; }
name=$1; shift
case $name in *[!A-Za-z0-9._-]*) echo "name: letters, digits, . _ - only"; exit 2 ;; esac
W=${PHYGAS_WORK:-$HOME/phygas-work}
STUB=$here/boinc_stub
OPT=${OPT:--O2}
CXXSTD=${CXXSTD--std=c++17}          # ${x-y}: an empty CXXSTD= stays empty
ALLOC=${ALLOC:-system}
if [ "$(uname)" = Darwin ]; then
  CC=${CC:-clang}; CXX=${CXX:-clang++}
  JEMALLOC=${JEMALLOC:-$(brew --prefix 2>/dev/null || echo /opt/homebrew)/lib/libjemalloc.a}
  GPULIBS="-lm -lpthread"
  BACKENDS=${BACKENDS:-cpu}
else
  CC=${CC:-cc}; CXX=${CXX:-c++}
  JEMALLOC=${JEMALLOC:--ljemalloc}
  GPULIBS="-lm -lpthread -ldl -lrt"
  if [ -z "${BACKENDS:-}" ]; then
    BACKENDS=cpu
    if command -v "${NVCC:-nvcc}" > /dev/null; then BACKENDS="cpu cuda"; fi
  fi
fi
case $ALLOC in
  system) L="" ;;
  jemalloc) L=$JEMALLOC ;;
  mimalloc) L=${MIMALLOC:--lmimalloc} ;;
  *) echo "ALLOC must be system, jemalloc or mimalloc"; exit 2 ;;
esac
for b in $BACKENDS; do
  case $b in cpu|cuda|opencl) ;; *) echo "BACKENDS: cpu, cuda or opencl (not $b)"; exit 2 ;; esac
done

T=$W/build/$name; out=$W/bins/$name
rm -rf "$T"; mkdir -p "$T/obj" "$out"
if [ -n "${REF:-}" ]; then
  git -C "$repo" rev-parse -q --verify "$REF^{commit}" > /dev/null || { echo "no git ref $REF in $repo"; exit 2; }
  git -C "$repo" archive "$REF" | tar -xf - -C "$T"
  srcdesc="$REF ($(git -C "$repo" rev-parse --short "$REF"))"
else
  gone=$(git -C "$repo" ls-files -d)
  [ -z "$gone" ] || { echo "tracked files are deleted in the working tree (restore them or git rm them):"; echo "$gone"; exit 2; }
  (cd "$repo" && git ls-files -z | tar --null -T - -cf -) | tar -xf - -C "$T"
  modified=$(git -C "$repo" status --porcelain --untracked-files=no | wc -l | tr -d ' ')
  srcdesc="working tree on $(git -C "$repo" rev-parse --short HEAD), $modified tracked files modified"
fi
[ -f "$T/src/main_phygas_boinc.cpp" ] && [ -f "$T/src/phygas.c" ] || { echo "$srcdesc has no src/main_phygas_boinc.cpp"; exit 1; }
for p in "$@"; do
  [ -f "$p" ] || { echo "no patch $p"; exit 2; }
  a="$(cd "$(dirname "$p")" && pwd)/$(basename "$p")"
  (cd "$T" && patch -p1 -s < "$a") || { echo "patch $p does not apply to $srcdesc"; exit 1; }
done

log=$T/build.log
: > "$log"
step() {   # one build command; its output goes to the log, a failure stops the script with the log's tail
  echo "+ $*" >> "$log"
  "$@" >> "$log" 2>&1 || { tail -n 25 "$log"; echo "build failed: $*"; echo "(whole log: $log)"; exit 1; }
}

# the C core is the same for every backend
CORE=""
for c in phygas coords potential montecarlo; do
  step $CC $OPT -std=c99 -I "$T/src" -c "$T/src/$c.c" -o "$T/obj/$c.o"
  CORE="$CORE $T/obj/$c.o"
done
step $CXX $OPT $CXXSTD -I "$STUB" -c "$STUB/boinc_stub.cpp" -o "$T/obj/stub.o"

for b in $BACKENDS; do
  rm -f "$out/$b"
  case $b in
    cpu)
      step $CXX $OPT $CXXSTD -I "$T/src" -I "$STUB" -c "$T/src/main_phygas_boinc.cpp" -o "$T/obj/main_cpu.o"
      step $CXX $OPT "$T/obj/main_cpu.o" "$T/obj/stub.o" $CORE $L -lm -lpthread -o "$out/cpu" ;;
    opencl)   # the OpenCL library is loaded at run time (cl_dynload.c), so there is no -lOpenCL
      step $CXX $OPT -std=c++17 -DPG_GPU -DPG_OPENCL -I "$T/src" -I "$STUB" -c "$T/src/main_phygas_boinc.cpp" -o "$T/obj/main_ocl.o"
      step $CC $OPT -std=c99 -DCL_TARGET_OPENCL_VERSION=120 -I "$T/src" -c "$T/src/cl_dynload.c" -o "$T/obj/cl_dynload.o"
      step $CC $OPT -std=c99 -D_POSIX_C_SOURCE=200809L -DCL_TARGET_OPENCL_VERSION=120 "-DEXPECTED_GPU_VENDOR=\"${VENDOR:-}\"" \
        -I "$T/src" -c "$T/src/phygas_opencl.c" -o "$T/obj/phygas_opencl.o"
      step $CXX $OPT "$T/obj/main_ocl.o" "$T/obj/stub.o" "$T/obj/cl_dynload.o" "$T/obj/phygas_opencl.o" $CORE $L $GPULIBS \
        -o "$out/opencl" ;;
    cuda)     # no nvcc, no CUDA app: this gap is in docs/BENCHMARK.md
      NVCC=${NVCC:-nvcc}
      command -v "$NVCC" > /dev/null || { echo "cuda: no $NVCC (CUDA toolkit); build the other backends or set NVCC"; exit 1; }
      CUDA_HOME=${CUDA_HOME:-$(cd "$(dirname "$(command -v "$NVCC")")/.." && pwd)}
      step "$NVCC" -O3 -std=c++17 --fmad=false -arch="${CUDA_ARCH:-native}" ${NVFLAGS:-} -I "$T/src" \
        -c "$T/src/phygas_cuda.cu" -o "$T/obj/phygas_cuda.o"
      step $CXX $OPT -std=c++17 -DPG_GPU -DPG_CUDA -I "$T/src" -I "$STUB" -I "$CUDA_HOME/include" \
        -c "$T/src/main_phygas_boinc.cpp" -o "$T/obj/main_cuda.o"
      step "$NVCC" -o "$out/cuda" "$T/obj/main_cuda.o" "$T/obj/stub.o" $CORE "$T/obj/phygas_cuda.o" ;;
  esac
  # MinGW-w64 appends .exe to an output name without one (cross_check.sh); keep the name the other scripts expect
  if [ -f "$out/$b.exe" ]; then mv -f "$out/$b.exe" "$out/$b"; fi
  [ -x "$out/$b" ] || { echo "build of $b left no executable at $out/$b, see $log"; exit 1; }
done

{
  echo "built $(date '+%Y-%m-%d %H:%M:%S') from $srcdesc"
  echo "backends: $BACKENDS   OPT=$OPT   CXXSTD=$CXXSTD   ALLOC=$ALLOC"
  echo "CC: $($CC --version 2>&1 | head -1)"
  echo "CXX: $($CXX --version 2>&1 | head -1)"
  if [ $# -gt 0 ]; then echo "patches: $*"; fi
} > "$out/build.info"
echo "built $out: $(cd "$out" && ls | grep -v build.info | tr '\n' ' ')($srcdesc; OPT=$OPT ALLOC=$ALLOC)"
