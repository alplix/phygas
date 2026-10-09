#ifndef NOSTOS_VEC3_H
#define NOSTOS_VEC3_H

#include <math.h>

typedef struct { double x, y, z; } Vec3;

// v3_make instead of (Vec3){...} compound literals: those are C-only, and this header is also
// included from C++ (the BOINC wrappers; MSVC's C++ rejects them). Same values, same arithmetic.
static inline Vec3 v3_make(double x, double y, double z) { Vec3 r; r.x = x; r.y = y; r.z = z; return r; }

static inline Vec3 v3_add(Vec3 a, Vec3 b) { return v3_make(a.x+b.x, a.y+b.y, a.z+b.z); }
static inline Vec3 v3_sub(Vec3 a, Vec3 b) { return v3_make(a.x-b.x, a.y-b.y, a.z-b.z); }
static inline Vec3 v3_scale(Vec3 a, double s) { return v3_make(a.x*s, a.y*s, a.z*s); }
static inline double v3_dot(Vec3 a, Vec3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static inline double v3_norm(Vec3 a) { return sqrt(v3_dot(a, a)); }

static inline Vec3 v3_cross(Vec3 a, Vec3 b) {
    return v3_make(
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x
    );
}

#endif
