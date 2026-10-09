// Phygas GPU backends (CUDA: phygas_cuda.cu, OpenCL: phygas_opencl.c) -- one interface, used by
// main_phygas_gpu_boinc.cpp. Device code: phygas_kernel_body.h (literal ports of the CPU screening).
//
// Every call is synchronous, so the caller can time a launch and size the next one (pg_next_slice):
// a launch must stay well under Windows' GPU watchdog (TDR, ~2 s) on weak-fp64 cards.
#ifndef PHYGAS_GPU_H
#define PHYGAS_GPU_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef void *PgBuf;

int pg_init(int platform, int device);   // 0 = ok
PgBuf pg_alloc(size_t bytes);            // NULL = failed
int pg_put(PgBuf b, const void *src, size_t bytes);
int pg_get(PgBuf b, void *dst, size_t bytes);
void pg_free(PgBuf b);

// Kernel launches (arguments as in phygas_kernel_body.h); 0 = ok.
int pg_lin(PgBuf ps0, PgBuf sig, PgBuf state, PgBuf trk, PgBuf params, double eps_pos, double eps_vel,
           double dt, int sub, int k0, int nk, int npts, int ns);
int pg_chi2(PgBuf trk, PgBuf csoa, PgBuf cn, PgBuf R2s, int nb, int nc, int npts, int p0, int np,
            PgBuf chi2, PgBuf kb, PgBuf sb);
int pg_finish(PgBuf trk, PgBuf csoa, PgBuf cn, PgBuf r50s, PgBuf rts, int nc, int npts, PgBuf hs, PgBuf hc,
              PgBuf hkb, PgBuf hsb, PgBuf hseed, int h0, int nh, PgBuf out);
int pg_pair(PgBuf ps, PgBuf steps, PgBuf dts, PgBuf st, PgBuf res, PgBuf params, int s0, int nk, int n);

#define PG_TARGET_LAUNCH_S 0.25

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
static inline double pg_now(void) {
    LARGE_INTEGER f, c;
    QueryPerformanceFrequency(&f);
    QueryPerformanceCounter(&c);
    return (double)c.QuadPart / (double)f.QuadPart;
}
#else
#include <time.h>
static inline double pg_now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + ts.tv_nsec * 1e-9;
}
#endif
// next launch size after one of size k took `elapsed` seconds: aim at PG_TARGET_LAUNCH_S, grow at
// most 4x per launch, stay within [1, max_k]
static inline int pg_next_slice(int k, double elapsed, int max_k) {
    double want = elapsed > 1e-4 ? k * (PG_TARGET_LAUNCH_S / elapsed) : k * 4.0;
    if (want > k * 4.0) want = k * 4.0;
    if (want > max_k) want = max_k;
    if (want < 1) want = 1;
    return (int)want;
}

#ifdef __cplusplus
}
#endif
#endif
