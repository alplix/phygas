#!/bin/bash
# Phygas GPU apps built on the bitboinc build server, from a directory containing src/:
#   Linux x86_64 CUDA      nvidia/cuda:12.8.1-devel-ubuntu20.04 container, static cudart
#   Linux x86_64 OpenCL    ubuntu:20.04 container, one binary per vendor (amd, nvidia)
#   Windows x86_64 OpenCL  host mingw-w64, static, per vendor
# (Windows CUDA is built with MSVC on the dev PC: scripts/build_gpu_windows.bat.)
# The Ubuntu 20.04 containers keep the GLIBC requirement at 2.31 (Debian 11 / Ubuntu 20.04 and newer); checked at the end.
# Intel GPUs are not targeted: most have no double precision, which Phygas needs.
set -e
V=${PHYGAS_VERSION:-2}
W=$PWD
mkdir -p out
BOINC_SRC=/opt/bitboinc-src/boinc
CSRC="coords potential montecarlo phygas"

cat > inner_boinc.sh <<'EOF'
export DEBIAN_FRONTEND=noninteractive
apt-get update -qq && apt-get install -y -qq build-essential m4 pkg-config file $EXTRA_PKGS >/dev/null
cp -r /boincsrc-orig /boincsrc && cd /boincsrc
for f in aclocal.m4 configure config.h.in Makefile.in Makefile; do find . -name "$f" -exec touch {} +; done
./configure --disable-client --disable-manager --disable-server --enable-apps --no-create --no-recursion --disable-maintainer-mode >/dev/null
make -C api -j4 >/dev/null && make -C lib -j4 >/dev/null
B=/boincsrc
cd /work
EOF

cat > inner_cuda.sh <<'EOF'
#!/bin/bash
set -e
. /work/inner_boinc.sh
G="-gencode arch=compute_50,code=sm_50 -gencode arch=compute_60,code=sm_60 -gencode arch=compute_70,code=sm_70 -gencode arch=compute_75,code=sm_75 -gencode arch=compute_80,code=sm_80 -gencode arch=compute_86,code=sm_86 -gencode arch=compute_89,code=sm_89 -gencode arch=compute_90,code=sm_90 -gencode arch=compute_100,code=sm_100 -gencode arch=compute_120,code=sm_120 -gencode arch=compute_120,code=compute_120"
mkdir -p obj_cuda
for c in coords potential montecarlo phygas; do gcc -O2 -std=c99 -I src -c src/$c.c -o obj_cuda/$c.o; done
nvcc -O3 -std=c++17 --fmad=false $G -I src -c src/phygas_cuda.cu -o obj_cuda/phygas_cuda.o
g++ -O2 -std=c++17 -DPG_GPU -DPG_CUDA -I src -I /usr/local/cuda/include -I $B/api -I $B/lib -c src/main_phygas_boinc.cpp -o obj_cuda/main.o
g++ obj_cuda/*.o $B/api/.libs/libboinc_api.a $B/lib/.libs/libboinc.a \
    -L /usr/local/cuda/targets/x86_64-linux/lib -lcudart_static -lculibos -lpthread -ldl -lrt -static-libstdc++ -static-libgcc \
    -o out/phygas_cuda_v${V}_linux_x86_64
EOF

cat > inner_opencl.sh <<'EOF'
#!/bin/bash
set -e
. /work/inner_boinc.sh
mkdir -p obj_ocl
for c in coords potential montecarlo phygas; do gcc -O2 -std=c99 -I src -c src/$c.c -o obj_ocl/$c.o; done
gcc -O2 -std=c99 -I src -DCL_TARGET_OPENCL_VERSION=120 -c src/cl_dynload.c -o obj_ocl/cl_dynload.o
g++ -O2 -std=c++17 -DPG_GPU -DPG_OPENCL -I src -I $B/api -I $B/lib -c src/main_phygas_boinc.cpp -o obj_ocl/main.o
for vt in ati:amd nvidia:nvidia; do
    tag=${vt%%:*}; vendor=${vt##*:}
    gcc -O2 -std=c99 -D_POSIX_C_SOURCE=200809L -I src -DCL_TARGET_OPENCL_VERSION=120 "-DEXPECTED_GPU_VENDOR=\"$vendor\"" -c src/phygas_opencl.c -o obj_ocl/phygas_opencl_$tag.o
    g++ obj_ocl/main.o obj_ocl/phygas_opencl_$tag.o obj_ocl/cl_dynload.o obj_ocl/coords.o obj_ocl/potential.o obj_ocl/montecarlo.o obj_ocl/phygas.o \
        $B/api/.libs/libboinc_api.a $B/lib/.libs/libboinc.a -lpthread -ldl -lrt -static-libstdc++ -static-libgcc \
        -o out/phygas_opencl_v${V}_linux_x86_64_$tag
done
EOF

echo "=== Linux CUDA ==="
sudo docker run --rm -e V=$V -e EXTRA_PKGS="" -v $BOINC_SRC:/boincsrc-orig:ro -v $W:/work nvidia/cuda:12.8.1-devel-ubuntu20.04 bash /work/inner_cuda.sh
echo "=== Linux OpenCL ==="
sudo docker run --rm -e V=$V -e EXTRA_PKGS="opencl-headers" -v $BOINC_SRC:/boincsrc-orig:ro -v $W:/work ubuntu:20.04 bash /work/inner_opencl.sh
sudo chown -R $(id -u):$(id -g) out obj_cuda obj_ocl 2>/dev/null || true

echo "=== Windows OpenCL (mingw) ==="
WB=/opt/bitboinc-src/boinc-mingw64
MG=/usr/bin/x86_64-w64-mingw32
mkdir -p obj_win
for c in $CSRC; do $MG-gcc -O2 -std=c99 -I src -c src/$c.c -o obj_win/$c.o; done
$MG-gcc -O2 -std=c99 -I src -DCL_TARGET_OPENCL_VERSION=120 -c src/cl_dynload.c -o obj_win/cl_dynload.o
$MG-g++ -O2 -std=c++17 -DPG_GPU -DPG_OPENCL -I src -I $WB/api -I $WB/lib -c src/main_phygas_boinc.cpp -o obj_win/main.o
for vt in ati:amd nvidia:nvidia; do
    tag=${vt%%:*}; vendor=${vt##*:}
    $MG-gcc -O2 -std=c99 -I src -DCL_TARGET_OPENCL_VERSION=120 "-DEXPECTED_GPU_VENDOR=\"$vendor\"" -c src/phygas_opencl.c -o obj_win/phygas_opencl_$tag.o
    $MG-g++ -static obj_win/main.o obj_win/phygas_opencl_$tag.o obj_win/cl_dynload.o obj_win/coords.o obj_win/potential.o obj_win/montecarlo.o obj_win/phygas.o \
        $WB/api/libboinc_api.a $WB/lib/libboinc.a -lpsapi -o out/phygas_opencl_v${V}_windows_x86_64_$tag.exe
done

for f in out/phygas_cuda_v${V}_linux_x86_64 out/phygas_opencl_v${V}_linux_x86_64_*; do
    echo "$f: needs GLIBC $(objdump -T $f | grep -oE 'GLIBC_[0-9.]+' | sort -V | tail -1); $(objdump -p $f | grep NEEDED | awk '{print $2}' | tr '\n' ' ')"
done
ls -la out/
