// Phygas CUDA backend (see phygas_gpu.h). Build with --fmad=false: the CPU builds do not contract
// a*b+c into FMAs, and keeping the same rounding keeps GPU and CPU results within rounding.
#include <cstdio>
#include <cuda_runtime.h>
#include "phygas_gpu.h"

#define PG_GLOBAL
#define PG_DEVICE __device__ static inline
#define PG_KERNEL __global__
#define PG_TID ((int)(blockIdx.x * blockDim.x + threadIdx.x))
#define PG_U64 unsigned long long
#define PG_LONG long long
#define PG_INF __longlong_as_double(0x7ff0000000000000LL)
#include "phygas_kernel_body.h"

#define CK(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { \
    fprintf(stderr, "phygas_cuda: %s failed: %s\n", #x, cudaGetErrorString(e_)); return 1; } } while (0)

static const int kThreads = 128;
static int blocks(int n) { return (n + kThreads - 1) / kThreads; }

extern "C" int pg_init(int, int device) {
    int n = 0;
    CK(cudaGetDeviceCount(&n));
    if (n <= 0) { fprintf(stderr, "phygas_cuda: no CUDA device\n"); return 1; }
    if (device < 0 || device >= n) device = 0;
    CK(cudaSetDevice(device));
    cudaDeviceProp prop;
    CK(cudaGetDeviceProperties(&prop, device));
    fprintf(stderr, "phygas_cuda: device %d: %s (sm_%d%d, %d SMs, %.0f MB)\n", device, prop.name, prop.major,
            prop.minor, prop.multiProcessorCount, prop.totalGlobalMem / 1048576.0);
    return 0;
}

extern "C" PgBuf pg_alloc(size_t bytes) {
    void *p = 0;
    if (cudaMalloc(&p, bytes ? bytes : 8) != cudaSuccess) return 0;
    return p;
}
extern "C" int pg_put(PgBuf b, const void *src, size_t bytes) { CK(cudaMemcpy(b, src, bytes, cudaMemcpyHostToDevice)); return 0; }
extern "C" int pg_get(PgBuf b, void *dst, size_t bytes) { CK(cudaMemcpy(dst, b, bytes, cudaMemcpyDeviceToHost)); return 0; }
extern "C" void pg_free(PgBuf b) { if (b) cudaFree(b); }

#define D(b) ((double *)(b))
#define I(b) ((int *)(b))

extern "C" int pg_lin(PgBuf ps0, PgBuf sig, PgBuf state, PgBuf trk, PgBuf params, double eps_pos, double eps_vel,
                      double dt, int sub, int k0, int nk, int npts, int ns) {
    pg_lin_k<<<blocks(ns), kThreads>>>(D(ps0), D(sig), D(state), D(trk), D(params), eps_pos, eps_vel, dt, sub, k0, nk, npts, ns);
    CK(cudaGetLastError());
    CK(cudaDeviceSynchronize());
    return 0;
}

extern "C" int pg_chi2(PgBuf trk, PgBuf csoa, PgBuf cn, PgBuf R2s, int nb, int nc, int npts, int p0, int np,
                       PgBuf chi2, PgBuf kb, PgBuf sb) {
    pg_chi2_k<<<blocks(np), kThreads>>>(D(trk), D(csoa), I(cn), D(R2s), nb, nc, npts, p0, np, D(chi2), I(kb), D(sb));
    CK(cudaGetLastError());
    CK(cudaDeviceSynchronize());
    return 0;
}

extern "C" int pg_finish(PgBuf trk, PgBuf csoa, PgBuf cn, PgBuf r50s, PgBuf rts, int nc, int npts, PgBuf hs, PgBuf hc,
                         PgBuf hkb, PgBuf hsb, PgBuf hseed, int h0, int nh, PgBuf out) {
    pg_finish_k<<<blocks(nh), kThreads>>>(D(trk), D(csoa), I(cn), D(r50s), D(rts), nc, npts, I(hs), I(hc), I(hkb),
                                          D(hsb), (const unsigned long long *)hseed, h0, nh, D(out));
    CK(cudaGetLastError());
    CK(cudaDeviceSynchronize());
    return 0;
}
