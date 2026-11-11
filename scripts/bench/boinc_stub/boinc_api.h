#ifndef PHYGAS_BENCH_BOINC_API_H
#define PHYGAS_BENCH_BOINC_API_H
// Minimal stand-in for BOINC's api/boinc_api.h, enough to build src/main_phygas_boinc.cpp without a BOINC checkout and
// run it as a plain process (scripts/bench/build_variant.sh, run_bench.sh). Only the calls the app makes are provided;
// boinc_stub.cpp says how each behaves. Written for C++98 on purpose (no default member initialisers): the CPU app's
// documented build line has no -std flag, and Apple clang 14 defaults to an old one.

struct BOINC_OPTIONS {
    int main_program, check_heartbeat, handle_process_control, send_status_msgs, direct_process_action;
    int multi_thread, multi_process, normal_thread_priority;
    BOINC_OPTIONS() { reset(); }
    void reset() {
        main_program = check_heartbeat = handle_process_control = send_status_msgs = direct_process_action = 0;
        multi_thread = multi_process = normal_thread_priority = 0;
    }
};

struct APP_INIT_DATA {
    int gpu_device_num;          // -1: not assigned, the app falls back to --device
    int gpu_opencl_dev_index;
    APP_INIT_DATA() : gpu_device_num(-1), gpu_opencl_dev_index(-1) {}
};

inline void boinc_options_defaults(BOINC_OPTIONS& o) { o.reset(); }
int boinc_init_options(BOINC_OPTIONS*);
int boinc_finish(int status);
int boinc_resolve_filename(const char* virtual_name, char* physical_name, int len);
int boinc_get_init_data(APP_INIT_DATA&);
int boinc_fraction_done(double);
int boinc_time_to_checkpoint();
int boinc_checkpoint_completed();
void boinc_begin_critical_section();
void boinc_end_critical_section();

#endif
