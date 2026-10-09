#!/bin/bash
# Phygas CPU BOINC app for every Linux and Windows platform bitboinc serves, on the bitboinc build
# server (cross compilers + per-target BOINC libs in /opt/bitboinc-src/boinc-<target>). Run from a
# directory containing src/ (this repo's). Mac and Android
# are built (as for Nostos) on their own hosts.
#
# Linux builds are fully static, so they carry no GLIBC version requirement at all (checked below:
# no GLIBC_ symbols). 32-bit x86 uses SSE2 arithmetic: with x87 80-bit intermediates results differ
# from every other platform (and Kampylos's fit even hung).
# Each target is independent: one failing doesn't stop the others.
V=${PHYGAS_VERSION:-1}
mkdir -p out
CSRC="phygas coords potential montecarlo"
LLVM=/opt/llvm-mingw/bin
X86_32="-msse2 -mfpmath=sse"

build() {  # tag  cc  cxx  boinc-dir  output  link-libs  [extra-flags]
    local tag=$1 cc="$2 $7" cxx="$3 $7" B=/opt/bitboinc-src/$4 outf=out/$5 libs=$6 o=obj_$1
    rm -rf $o; mkdir -p $o
    echo "=== $tag ==="
    local ok=1
    for c in $CSRC; do
        $cc -O2 -std=c99 -I src -c src/$c.c -o $o/$c.o || ok=0
    done
    [ $ok = 1 ] && $cxx -O2 -I src -I $B/api -I $B/lib -c src/main_phygas_boinc.cpp -o $o/main.o || ok=0
    [ $ok = 1 ] && $cxx -static $o/*.o $B/api/libboinc_api.a $B/lib/libboinc.a $libs -o $outf || ok=0
    if [ $ok = 1 ]; then
        file $outf | cut -c1-140
        case $outf in *.exe) ;; *) echo "  GLIBC-versioned symbols: $(objdump -T $outf 2>/dev/null | grep -c GLIBC_)";; esac
    else
        echo "$tag FAILED"
    fi
}

build x86_64  gcc                       g++                         boinc          phygas_v${V}_linux_x86_64     "-lpthread -lm"
build aarch64 aarch64-linux-gnu-gcc     aarch64-linux-gnu-g++       boinc-aarch64  phygas_v${V}_linux_aarch64    "-lpthread -lm"
build armhf   arm-linux-gnueabihf-gcc   arm-linux-gnueabihf-g++     boinc-armhf    phygas_v${V}_linux_armhf      "-lpthread -lm"
build i686    i686-linux-gnu-gcc        i686-linux-gnu-g++          boinc-i686     phygas_v${V}_linux_i686       "-lpthread -lm" "$X86_32"
build ppc64le powerpc64le-linux-gnu-gcc powerpc64le-linux-gnu-g++   boinc-ppc64le  phygas_v${V}_linux_ppc64le    "-lpthread -lm"
build riscv64 riscv64-linux-gnu-gcc     riscv64-linux-gnu-g++       boinc-riscv64  phygas_v${V}_linux_riscv64    "-lpthread -lm"
build win64   /usr/bin/x86_64-w64-mingw32-gcc /usr/bin/x86_64-w64-mingw32-g++ boinc-mingw64 phygas_v${V}_windows_x86_64.exe "-lpsapi"
build win32   /usr/bin/i686-w64-mingw32-gcc   /usr/bin/i686-w64-mingw32-g++   boinc-mingw32 phygas_v${V}_windows_x86.exe    "-lpsapi" "$X86_32"
build winarm64 $LLVM/aarch64-w64-mingw32-clang $LLVM/aarch64-w64-mingw32-clang++ boinc-winarm64 phygas_v${V}_windows_arm64.exe "-lpsapi"
ls -la out/
