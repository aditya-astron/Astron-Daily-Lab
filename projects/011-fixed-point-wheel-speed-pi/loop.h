/*
 * Astron Built this AP
 * Project 011 | Embedded firmware | Fixed-point PI speed loop for a reaction wheel with anti-windup
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-03 IST
 */
/* Closed-loop harness: wheel parameters, pole-placement gain design and a
 * sampled-data simulation of the Q16.16 controller driving the wheel. */
#ifndef LOOP_H
#define LOOP_H

#include "pi_q16.h"
#include "wheel.h"

typedef struct {
    double J, kt, b;   /* plant */
    double i_max;      /* motor driver current limit, A */
    double ts;         /* controller period, s */
    double wn, zeta;   /* closed-loop design targets */
} loop_cfg_t;

typedef struct { double kp, ki; } pi_gains_t;

loop_cfg_t loop_default_cfg(void);
pi_gains_t loop_design(const loop_cfg_t *cfg);

/* Analytic unsaturated step response of the continuous loop, per unit step. */
double loop_analytic_step(const loop_cfg_t *cfg, const pi_gains_t *g, double t);

typedef double (*profile_fn)(double t);

/* Runs n samples. ref(t), load(t) give speed reference (rad/s) and load
 * torque (N m). Writes speed and commanded current per sample (arrays of n,
 * may be NULL). */
void loop_run(const loop_cfg_t *cfg, const pi_gains_t *g, aw_mode_t aw,
              double w0, profile_fn ref, profile_fn load, int n,
              double *w_out, double *i_out);

/* Same controller law in double precision, used to bound the fixed-point
 * quantization error. Mirrors pi_q16_step exactly apart from arithmetic. */
typedef struct { double kp, ki_dt, u_max, integ; aw_mode_t aw; } pi_dbl_t;
void pi_dbl_init(pi_dbl_t *c, double kp, double ki, double ts, double u_max,
                 aw_mode_t aw);
double pi_dbl_step(pi_dbl_t *c, double ref, double meas);

#endif
