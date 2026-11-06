/*
 * Astron Built this AP
 * Project 017 | Attitude control | Rate- and torque-limited eigenaxis slew with reaction wheels
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-06 IST
 *
 * Shared reference scenario: a ~1.2 t agile imaging satellite with a slightly
 * non-diagonal inertia, three body-axis wheels, and bias momentum already
 * stored in the wheels from earlier disturbance rejection.
 */
#ifndef SCENARIO_H
#define SCENARIO_H

#include "slew.h"

static inline spacecraft scenario_spacecraft(void)
{
    spacecraft sc = {
        .J = {{{1200.0, -15.0, 22.0},
               {-15.0, 1500.0, -8.0},
               {22.0, -8.0, 900.0}}},
        .tau_max = 0.6,
        .h_max = 16.0,
    };
    m3_inverse(&sc.J, &sc.Jinv);
    return sc;
}

static inline vec3 scenario_h0(void) { return v3(1.5, -2.0, 0.8); }

#define SCENARIO_ZETA   0.9
#define SCENARIO_WN     0.1   /* rad/s, well below the first panel mode */
#define SCENARIO_MARGIN 0.8
#define SCENARIO_BRAKE  0.6

#endif
