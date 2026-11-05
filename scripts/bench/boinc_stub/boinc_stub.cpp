// Fake BOINC runtime for running the Phygas app as a plain process (see boinc_api.h). Files are plain files in the
// working directory ("in", "clusters", "out", "phygas_checkpoint"). boinc_time_to_checkpoint() follows the real API's
// contract (1 at most once per period, re-armed by boinc_checkpoint_completed(); each 1 opens a critical section);
// the period comes from PHYGAS_STUB_CHECKPOINT_SECONDS (default 60, the client's default; 0 = at every call).
//
// Left out on purpose, so their cost is not in the benchmark numbers: the timer thread and its SIGALRM, heartbeat,
// status messages through shared memory, and the 2 s sleep in the real boinc_finish(). boinc_fraction_done() keeps
// the value in memory (the real one only sets a variable too) and boinc_finish() writes it to stub_fraction_done,
// so that a progress report does not cost a file write per target.
#include "boinc_api.h"
#include "filesys.h"
#include <cstdio>
#include <cstdlib>
#include <time.h>

static double stub_period = 60.0;
static double stub_last = 0.0;
static double stub_fraction = 0.0;
static int stub_critical = 0;
static long stub_checkpoints = 0;

static double stub_now() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + ts.tv_nsec * 1e-9;
}

int boinc_init_options(BOINC_OPTIONS*) {
    const char* p = getenv("PHYGAS_STUB_CHECKPOINT_SECONDS");
    if (p) stub_period = atof(p);
    stub_last = stub_now();
    return 0;
}

int boinc_finish(int status) {
    FILE* fp = fopen("stub_fraction_done", "w");
    if (fp) { fprintf(fp, "%.6f\n", stub_fraction); fclose(fp); }
    fprintf(stderr, "boinc_stub: boinc_finish(%d), %ld checkpoints reported\n", status, stub_checkpoints);
    fflush(stderr);
    exit(status);
}

int boinc_resolve_filename(const char* virtual_name, char* physical_name, int len) {
    snprintf(physical_name, len, "%s", virtual_name);
    return 0;
}

// no GPU assignment: the GPU apps take --device (and --cl-platform), as they do when a task carries no init_data
int boinc_get_init_data(APP_INIT_DATA& d) {
    d = APP_INIT_DATA();
    return 0;
}

int boinc_fraction_done(double f) {
    stub_fraction = f;
    return 0;
}

int boinc_time_to_checkpoint() {
    if (stub_now() - stub_last < stub_period) return 0;
    boinc_begin_critical_section();
    return 1;
}

int boinc_checkpoint_completed() {
    stub_last = stub_now();
    stub_checkpoints++;
    boinc_end_critical_section();
    return 0;
}

void boinc_begin_critical_section() { stub_critical++; }
void boinc_end_critical_section() { if (stub_critical > 0) stub_critical--; }

FILE* boinc_fopen(const char* path, const char* mode) { return fopen(path, mode); }
