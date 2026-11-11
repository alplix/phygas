#ifndef PHYGAS_BENCH_FILESYS_H
#define PHYGAS_BENCH_FILESYS_H
// Stand-in for BOINC's lib/filesys.h (see boinc_api.h in this directory).
#include <cstdio>

FILE* boinc_fopen(const char* path, const char* mode);

#endif
