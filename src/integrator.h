// Orbit integration: leapfrog (kick-drift-kick), symplectic, standard choice
// for orbit integration in a smooth background potential -- exact energy and
// angular momentum conservation is impossible for any finite-step method, but
// leapfrog's error oscillates around the true value rather than drifting, so
// long integrations stay bounded instead of secularly decaying. See
// tests/test_kepler.c for the conservation + closed-orbit checks this
// property is validated against.
#ifndef NOSTOS_INTEGRATOR_H
#define NOSTOS_INTEGRATOR_H

#include "potential.h"
#include "vec3.h"

typedef struct {
    Vec3 pos; // kpc
    Vec3 vel; // kpc/Myr
} State;

// Advances `state` by exactly `dt` (Myr; negative dt integrates backward --
// leapfrog is time-symmetric, so this works without any special-casing).
static inline void leapfrog_step(const Potential *pot, State *state, double dt) {
    Vec3 a0 = potential_acceleration(pot, state->pos);
    // kick (half step)
    Vec3 v_half = v3_add(state->vel, v3_scale(a0, dt * 0.5));
    // drift (full step)
    state->pos = v3_add(state->pos, v3_scale(v_half, dt));
    // kick (half step), using the acceleration at the NEW position
    Vec3 a1 = potential_acceleration(pot, state->pos);
    state->vel = v3_add(v_half, v3_scale(a1, dt * 0.5));
}

// Integrates from state's current value over `total_time` (Myr, may be
// negative for backward integration) using `n_steps` fixed-size leapfrog
// steps. Calls `callback` (if non-NULL) after every step with the current
// time-since-start and state, e.g. to record a trajectory or check for a
// closest approach.
typedef void (*IntegratorCallback)(double t, const State *state, void *user_data);

static inline void integrate(
    const Potential *pot, State *state, double total_time, int n_steps,
    IntegratorCallback callback, void *user_data
) {
    double dt = total_time / n_steps;
    for (int i = 0; i < n_steps; i++) {
        leapfrog_step(pot, state, dt);
        if (callback) callback((i + 1) * dt, state, user_data);
    }
}

#endif
