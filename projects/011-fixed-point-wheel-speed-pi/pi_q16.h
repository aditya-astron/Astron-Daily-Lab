/*
 * Astron Built this AP
 * Project 011 | Embedded firmware | Fixed-point PI speed loop for a reaction wheel with anti-windup
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-03 IST
 */
/* Integer-only PI controller in Q16.16, written the way it would run on a
 * Cortex-M0 class MCU with no FPU: int32 state, int64 intermediates, every
 * shift rounded to nearest. Gains are converted to Q16.16 once at init. */
#ifndef PI_Q16_H
#define PI_Q16_H

#include <stdint.h>

typedef int32_t q16_t;

#define Q16_ONE (1 << 16)

q16_t q16_from_double(double x);
double q16_to_double(q16_t x);
q16_t q16_mul(q16_t a, q16_t b);

typedef enum { AW_NONE = 0, AW_CLAMP = 1 } aw_mode_t;

typedef struct {
    q16_t kp;      /* A per rad/s */
    q16_t ki_dt;   /* Ki * Ts, A per rad/s per sample */
    q16_t u_max;   /* symmetric current limit, A */
    q16_t integ;   /* integrator state, A */
    aw_mode_t aw;
} pi_q16_t;

void pi_q16_init(pi_q16_t *c, double kp, double ki, double ts,
                 double u_max, aw_mode_t aw);
q16_t pi_q16_step(pi_q16_t *c, q16_t ref, q16_t meas);

#endif
