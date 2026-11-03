/*
 * Astron Built this AP
 * Project 011 | Embedded firmware | Fixed-point PI speed loop for a reaction wheel with anti-windup
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-03 IST
 */
/* Verification of the fixed-point speed loop:
 *  - Q16.16 multiply rounding, including negative operands
 *  - gains place the closed-loop poles where designed
 *  - small-signal response matches the analytic continuous response
 *  - integer controller tracks a double-precision twin to a few LSB
 *  - saturated slew time matches the closed-form current-limited ramp
 *  - anti-windup removes the windup overshoot
 *  - load-step dip matches the analytic disturbance response, zero
 *    steady-state error */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "loop.h"

static int fails = 0;

static void check(int ok, const char *what, double got, double want)
{
    printf("[%s] %-46s got %-12.6g ref %.6g\n", ok ? "PASS" : "FAIL", what,
           got, want);
    if (!ok) fails++;
}

static double ref_101(double t) { (void)t; return 101.0; }
static double ref_300(double t) { (void)t; return 300.0; }
static double load_5m(double t) { return t >= 0.5 ? 5.0e-3 : 0.0; }

int main(void)
{
    loop_cfg_t c = loop_default_cfg();
    pi_gains_t g = loop_design(&c);

    /* Q16.16 arithmetic */
    q16_t prod = q16_mul(q16_from_double(1.5), q16_from_double(-2.25));
    check(prod == q16_from_double(-3.375), "q16_mul 1.5 * -2.25",
          q16_to_double(prod), -3.375);
    check(q16_mul(1, Q16_ONE / 2) == 1, "q16_mul rounds 0.5 LSB up",
          q16_mul(1, Q16_ONE / 2), 1);
    check(q16_mul(-3, Q16_ONE / 2) == -1, "q16_mul -1.5 LSB -> -1",
          q16_mul(-3, Q16_ONE / 2), -1);

    /* Pole placement: J s^2 + (b + Kt Kp) s + Kt Ki */
    double a1 = (c.b + c.kt * g.kp) / c.J, a0 = c.kt * g.ki / c.J;
    check(fabs(a1 - 2 * c.zeta * c.wn) < 1e-9, "char. poly s^1 coeff = 2 zeta wn",
          a1, 2 * c.zeta * c.wn);
    check(fabs(a0 - c.wn * c.wn) < 1e-9, "char. poly s^0 coeff = wn^2", a0,
          c.wn * c.wn);

    /* Small-signal step vs continuous analytic. The sampled loop differs
     * from the continuous one by O(wn Ts); allow 3 % of the step. */
    int n = 400;
    double *w = malloc(n * sizeof *w), *iq = malloc(n * sizeof *iq);
    loop_run(&c, &g, AW_CLAMP, 100.0, ref_101, NULL, n, w, iq);
    double err = 0, imax_seen = 0;
    for (int k = 0; k < n; k++) {
        err = fmax(err, fabs((w[k] - 100.0) - loop_analytic_step(&c, &g, k * c.ts)));
        imax_seen = fmax(imax_seen, fabs(iq[k]));
    }
    check(imax_seen < c.i_max, "small step stays unsaturated", imax_seen, c.i_max);
    check(err < 0.03, "small step max |sim - analytic| < 0.03", err, 0.03);
    check(fabs(w[n - 1] - 101.0) < 1e-3, "small step final value", w[n - 1], 101.0);

    /* Fixed-point vs double controller fed the identical measurement
     * sequence (the recorded closed-loop speeds). */
    pi_dbl_t d;
    pi_dbl_init(&d, g.kp, g.ki, c.ts, c.i_max, AW_CLAMP);
    d.integ = q16_to_double(q16_from_double(c.b * 100.0 / c.kt));
    double du = 0;
    for (int k = 0; k < n; k++) {
        double meas = q16_to_double(q16_from_double(w[k]));
        du = fmax(du, fabs(pi_dbl_step(&d, 101.0, meas) - iq[k]));
    }
    check(du < 20.0 / Q16_ONE, "Q16.16 vs double output < 20 LSB", du * Q16_ONE, 20);
    free(w); free(iq);

    /* Large step: current-limited ramp then capture. */
    n = 10000;
    w = malloc(n * sizeof *w);
    double w_inf = c.kt * c.i_max / c.b;
    double t90_a = c.J / c.b * log(w_inf / (w_inf - 270.0));
    double peak[2];
    for (int m = 0; m < 2; m++) {
        aw_mode_t aw = m ? AW_CLAMP : AW_NONE;
        loop_run(&c, &g, aw, 0.0, ref_300, NULL, n, w, NULL);
        double t90 = -1;
        peak[m] = 0;
        for (int k = 0; k < n; k++) {
            if (t90 < 0 && w[k] >= 270.0) t90 = k * c.ts;
            peak[m] = fmax(peak[m], w[k]);
        }
        if (m) check(fabs(t90 - t90_a) <= 2 * c.ts, "saturated t90 vs analytic ramp",
                     t90, t90_a);
        if (m) check(fabs(w[n - 1] - 300.0) < 1e-3, "large step final value",
                     w[n - 1], 300.0);
    }
    check((peak[0] - 300.0) / 300.0 > 0.5, "no anti-windup overshoot > 50 %",
          100 * (peak[0] - 300.0) / 300.0, 50);
    check((peak[1] - 300.0) / 300.0 < 0.01, "clamp anti-windup overshoot < 1 %",
          100 * (peak[1] - 300.0) / 300.0, 1);
    free(w);

    /* Load torque step: W/Td = -s / (J(s^2 + 2 zeta wn s + wn^2)). */
    n = 2500;
    w = malloc(n * sizeof *w);
    loop_run(&c, &g, AW_CLAMP, 300.0, ref_300, load_5m, n, w, NULL);
    double dip = 0;
    for (int k = 0; k < n; k++) dip = fmin(dip, w[k] - 300.0);
    double sg = c.zeta * c.wn, wd = c.wn * sqrt(1 - c.zeta * c.zeta);
    double tp = atan(wd / sg) / wd;
    double dip_a = -5.0e-3 / (c.J * wd) * exp(-sg * tp) * sin(wd * tp);
    check(fabs(dip / dip_a - 1) < 0.05, "load-step dip within 5 % of analytic",
          dip, dip_a);
    check(fabs(w[n - 1] - 300.0) < 1e-3, "zero steady-state error under load",
          w[n - 1] - 300.0, 0);
    free(w);

    printf("\n%s (%d failure%s)\n", fails ? "FAILED" : "ALL PASSED", fails,
           fails == 1 ? "" : "s");
    return fails != 0;
}
