/*
 * Astron Built this AP
 * Project 055 | Embedded firmware | Fixed-point CORDIC for a coarse sun sensor
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-21 IST
 *
 * Integer-only CORDIC for a 32-bit MCU without an FPU.
 *
 * Angles are binary angle measurement (BAM): a uint32_t where 2^32 is one full
 * turn, so wraparound is free and signed reinterpretation gives (-pi, pi].
 * Iteration count N (1..CORDIC_MAX_ITER) is a run-time argument so the same
 * code can trade cycles for accuracy.
 *
 * Assumes arithmetic right shift of negative int32_t (true on GCC/Clang for
 * every ARM and RISC-V target; checked in the test).
 */
#ifndef CORDIC_H
#define CORDIC_H

#include <stdint.h>

#define CORDIC_MAX_ITER 30
#define BAM_PER_TURN 4294967296.0          /* 2^32, for host-side conversions */

extern const int32_t cordic_atan_bam[CORDIC_MAX_ITER];  /* atan(2^-i) in BAM */
extern const int32_t cordic_inv_gain_q31[CORDIC_MAX_ITER]; /* 1/K_N, N = idx+1 */

/* Rotation mode: cos and sin of a BAM angle, Q2.30 (1.0 = 2^30). */
void cordic_sincos(uint32_t angle, int n_iter, int32_t *cos_q30, int32_t *sin_q30);

/*
 * Vectoring mode: angle = atan2(y, x) in BAM and mag = sqrt(x^2 + y^2) in the
 * same units as the inputs, rounded. Any int32_t inputs are accepted,
 * including INT32_MIN; inputs are block-normalised so the larger component
 * sits in [2^28, 2^29) before iterating, which keeps the datapath error fixed
 * relative to the vector length whatever the input scale.
 */
void cordic_atan2_mag(int32_t x, int32_t y, int n_iter, uint32_t *angle, uint32_t *mag);

#endif
