// Phygas OpenCL backend (see phygas_gpu.h). Same device code as the CUDA backend
// (phygas_kernel_body.h, embedded by scripts/embed_kernel.py), double precision required.
// The OpenCL library is loaded at run time (cl_dynload.c, shared with bitboinc's apps); per-vendor
// binaries are built with -DEXPECTED_GPU_VENDOR="amd"/"nvidia" so a task never runs on the wrong GPU.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cl_dynload.h"
#include "phygas_gpu.h"
#include "phygas_kernel_cl.inc"

#ifndef EXPECTED_GPU_VENDOR
#define EXPECTED_GPU_VENDOR ""
#endif

static cl_context g_ctx;
static cl_command_queue g_q;
static cl_program g_prog;
static cl_kernel g_lin, g_chi2, g_fin;
static cl_device_id g_dev;
static size_t g_local = 64;

static cl_kernel make_kernel(const char *name) {
    cl_int err;
    cl_kernel k = cl_api.CreateKernel(g_prog, name, &err);
    if (!k) { fprintf(stderr, "phygas_opencl: clCreateKernel(%s) failed (%d)\n", name, err); return NULL; }
    size_t wg = 0;
    cl_api.GetKernelWorkGroupInfo(k, g_dev, CL_KERNEL_WORK_GROUP_SIZE, sizeof wg, &wg, NULL);
    if (wg && wg < g_local) g_local = wg;
    return k;
}

int pg_init(int platform_id, int device_id) {
    if (!cl_dynload_init()) { fprintf(stderr, "phygas_opencl: could not load an OpenCL library\n"); return 1; }
    cl_uint np = 0;
    if (cl_api.GetPlatformIDs(0, NULL, &np) != CL_SUCCESS || np == 0) { fprintf(stderr, "phygas_opencl: no OpenCL platforms\n"); return 1; }
    cl_platform_id *ps = (cl_platform_id *)malloc(np * sizeof *ps);
    cl_api.GetPlatformIDs(np, ps, NULL);
    // BOINC's --device is an index within the GPU's platform, which need not be platform 0: take the
    // given platform if it matches the vendor and has devices, otherwise the first one that does.
    cl_platform_id plat = 0;
    cl_uint nd = 0;
    if (platform_id >= 0 && (cl_uint)platform_id < np && cl_platform_matches_vendor(ps[platform_id], EXPECTED_GPU_VENDOR)) {
        cl_api.GetDeviceIDs(ps[platform_id], CL_DEVICE_TYPE_GPU, 0, NULL, &nd);
        if (nd > 0) plat = ps[platform_id];
    }
    for (cl_uint i = 0; !plat && i < np; i++) {
        if (!cl_platform_matches_vendor(ps[i], EXPECTED_GPU_VENDOR)) continue;
        cl_uint n = 0;
        if (cl_api.GetDeviceIDs(ps[i], CL_DEVICE_TYPE_GPU, 0, NULL, &n) == CL_SUCCESS && n > 0) { plat = ps[i]; nd = n; }
    }
    free(ps);
    if (!plat) { fprintf(stderr, "phygas_opencl: no \"%s\" OpenCL GPU found\n", EXPECTED_GPU_VENDOR); return 1; }
    cl_device_id *ds = (cl_device_id *)malloc(nd * sizeof *ds);
    cl_api.GetDeviceIDs(plat, CL_DEVICE_TYPE_GPU, nd, ds, NULL);
    if (device_id < 0 || (cl_uint)device_id >= nd) device_id = 0;
    device_id = cl_pick_vendor_device(ds, nd, device_id, EXPECTED_GPU_VENDOR);
    g_dev = ds[device_id];
    free(ds);
    char name[256] = {0};
    cl_api.GetDeviceInfo(g_dev, CL_DEVICE_NAME, sizeof name, name, NULL);
    cl_device_fp_config fp64 = 0;
    cl_api.GetDeviceInfo(g_dev, CL_DEVICE_DOUBLE_FP_CONFIG, sizeof fp64, &fp64, NULL);
    cl_ulong mem = 0;
    cl_api.GetDeviceInfo(g_dev, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof mem, &mem, NULL);
    fprintf(stderr, "phygas_opencl: using \"%s\" (%.0f MB)\n", name, mem / 1048576.0);
    if (!fp64) { fprintf(stderr, "phygas_opencl: this GPU has no double precision (fp64), which Phygas needs\n"); return 1; }

    cl_int err;
    g_ctx = cl_api.CreateContext(NULL, 1, &g_dev, NULL, NULL, &err);
    if (!g_ctx) { fprintf(stderr, "phygas_opencl: clCreateContext failed (%d)\n", err); return 1; }
    g_q = cl_api.CreateCommandQueue(g_ctx, g_dev, 0, &err);
    if (!g_q) { fprintf(stderr, "phygas_opencl: clCreateCommandQueue failed (%d)\n", err); return 1; }
    const char *src = kPhygasKernelCL;
    size_t len = strlen(src);
    g_prog = cl_api.CreateProgramWithSource(g_ctx, 1, &src, &len, &err);
    err = g_prog ? cl_api.BuildProgram(g_prog, 1, &g_dev, "-cl-std=CL1.2", NULL, NULL) : err;
    if (err != CL_SUCCESS) {
        size_t n = 0;
        cl_api.GetProgramBuildInfo(g_prog, g_dev, CL_PROGRAM_BUILD_LOG, 0, NULL, &n);
        char *log = (char *)malloc(n + 1);
        if (log) {
            cl_api.GetProgramBuildInfo(g_prog, g_dev, CL_PROGRAM_BUILD_LOG, n, log, NULL);
            log[n] = 0;
            fprintf(stderr, "phygas_opencl: kernel build failed (%d):\n%s\n", err, log);
            free(log);
        }
        return 1;
    }
    g_lin = make_kernel("pg_lin_k");
    g_chi2 = make_kernel("pg_chi2_k");
    g_fin = make_kernel("pg_finish_k");
    return g_lin && g_chi2 && g_fin ? 0 : 1;
}

PgBuf pg_alloc(size_t bytes) {
    cl_int err;
    cl_mem m = cl_api.CreateBuffer(g_ctx, CL_MEM_READ_WRITE, bytes ? bytes : 8, NULL, &err);
    return err == CL_SUCCESS ? (PgBuf)m : NULL;
}
int pg_put(PgBuf b, const void *src, size_t bytes) {
    if (!bytes) return 0;
    cl_int err = cl_api.EnqueueWriteBuffer(g_q, (cl_mem)b, CL_TRUE, 0, bytes, src, 0, NULL, NULL);
    if (err != CL_SUCCESS) { fprintf(stderr, "phygas_opencl: write failed (%d)\n", err); return 1; }
    return 0;
}
int pg_get(PgBuf b, void *dst, size_t bytes) {
    if (!bytes) return 0;
    cl_int err = cl_api.EnqueueReadBuffer(g_q, (cl_mem)b, CL_TRUE, 0, bytes, dst, 0, NULL, NULL);
    if (err != CL_SUCCESS) { fprintf(stderr, "phygas_opencl: read failed (%d)\n", err); return 1; }
    return 0;
}
void pg_free(PgBuf b) { if (b) cl_api.ReleaseMemObject((cl_mem)b); }

// argument helpers: buffers, doubles and ints in kernel-signature order
#define ARG_B(k, i, b) do { cl_mem m_ = (cl_mem)(b); cl_api.SetKernelArg(k, i, sizeof(cl_mem), &m_); } while (0)
#define ARG_D(k, i, v) do { double d_ = (v); cl_api.SetKernelArg(k, i, sizeof(double), &d_); } while (0)
#define ARG_I(k, i, v) do { cl_int i_ = (v); cl_api.SetKernelArg(k, i, sizeof(cl_int), &i_); } while (0)

static int run(cl_kernel k, int n, const char *what) {
    if (n <= 0) return 0;
    size_t global = ((size_t)n + g_local - 1) / g_local * g_local;
    cl_int err = cl_api.EnqueueNDRangeKernel(g_q, k, 1, NULL, &global, &g_local, 0, NULL, NULL);
    if (err == CL_SUCCESS) err = cl_api.Finish(g_q);
    if (err != CL_SUCCESS) { fprintf(stderr, "phygas_opencl: %s kernel failed (%d)\n", what, err); return 1; }
    return 0;
}

int pg_lin(PgBuf ps0, PgBuf sig, PgBuf state, PgBuf trk, PgBuf params, double eps_pos, double eps_vel,
           double dt, int sub, int k0, int nk, int npts, int ns) {
    cl_kernel k = g_lin;
    ARG_B(k, 0, ps0); ARG_B(k, 1, sig); ARG_B(k, 2, state); ARG_B(k, 3, trk); ARG_B(k, 4, params);
    ARG_D(k, 5, eps_pos); ARG_D(k, 6, eps_vel); ARG_D(k, 7, dt);
    ARG_I(k, 8, sub); ARG_I(k, 9, k0); ARG_I(k, 10, nk); ARG_I(k, 11, npts); ARG_I(k, 12, ns);
    return run(k, ns, "lin");
}

int pg_chi2(PgBuf trk, PgBuf csoa, PgBuf cn, PgBuf R2s, int nb, int nc, int npts, int p0, int np,
            PgBuf chi2, PgBuf kb, PgBuf sb) {
    cl_kernel k = g_chi2;
    ARG_B(k, 0, trk); ARG_B(k, 1, csoa); ARG_B(k, 2, cn); ARG_B(k, 3, R2s);
    ARG_I(k, 4, nb); ARG_I(k, 5, nc); ARG_I(k, 6, npts); ARG_I(k, 7, p0); ARG_I(k, 8, np);
    ARG_B(k, 9, chi2); ARG_B(k, 10, kb); ARG_B(k, 11, sb);
    return run(k, np, "chi2");
}

int pg_finish(PgBuf trk, PgBuf csoa, PgBuf cn, PgBuf r50s, PgBuf rts, int nc, int npts, PgBuf hs, PgBuf hc,
              PgBuf hkb, PgBuf hsb, PgBuf hseed, int h0, int nh, PgBuf out) {
    cl_kernel k = g_fin;
    ARG_B(k, 0, trk); ARG_B(k, 1, csoa); ARG_B(k, 2, cn); ARG_B(k, 3, r50s); ARG_B(k, 4, rts);
    ARG_I(k, 5, nc); ARG_I(k, 6, npts);
    ARG_B(k, 7, hs); ARG_B(k, 8, hc); ARG_B(k, 9, hkb); ARG_B(k, 10, hsb); ARG_B(k, 11, hseed);
    ARG_I(k, 12, h0); ARG_I(k, 13, nh); ARG_B(k, 14, out);
    return run(k, nh, "finish");
}
