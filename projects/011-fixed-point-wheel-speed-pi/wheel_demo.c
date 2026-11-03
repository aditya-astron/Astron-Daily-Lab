/*
 * Astron Built this AP
 * Project 011 | Embedded firmware | Fixed-point PI speed loop for a reaction wheel with anti-windup
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-03 IST
 */
/* Runs the three scenarios that size a wheel speed loop: small-signal
 * tracking, a large saturated speed change with and without anti-windup,
 * and rejection of a load torque step. Prints metrics against the
 * continuous-time predictions. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "loop.h"

#define W_BIAS 100.0
#define DW_SMALL 1.0
#define W_BIG 300.0
#define TAU_STEP 5.0e-3
#define T_LOAD 0.5

static double ref_small(double t) { (void)t; return W_BIAS + DW_SMALL; }
static double ref_big(double t) { (void)t; return W_BIG; }
static double ref_hold(double t) { (void)t; return W_BIG; }
static double load_step(double t) { return t >= T_LOAD ? TAU_STEP : 0.0; }

static void big_step_metrics(const loop_cfg_t *c, const pi_gains_t *g,
                             aw_mode_t aw, const char *name)
{
    int n = (int)(10.0 / c->ts);
    double *w = malloc(n * sizeof *w);
    loop_run(c, g, aw, 0.0, ref_big, NULL, n, w, NULL);

    double peak = 0, t90 = -1, tset = 0;
    for (int k = 0; k < n; k++) {
        if (w[k] > peak) peak = w[k];
        if (t90 < 0 && w[k] >= 0.9 * W_BIG) t90 = k * c->ts;
        if (fabs(w[k] - W_BIG) > 0.01 * W_BIG) tset = (k + 1) * c->ts;
    }
    printf("   %-8s overshoot %6.2f %%   t90 %.3f s   ", name,
           100.0 * (peak - W_BIG) / W_BIG, t90);
    if (tset >= n * c->ts)
        printf("not settled in %.0f s (saturated limit cycle)\n", n * c->ts);
    else
        printf("1%% settle %.3f s\n", tset);
    free(w);
}

int main(void)
{
    loop_cfg_t c = loop_default_cfg();
    pi_gains_t g = loop_design(&c);
    printf("Wheel: J=%.1e kg m^2  Kt=%.3f N m/A  b=%.1e N m s  Imax=%.2f A  Ts=%.0f us\n",
           c.J, c.kt, c.b, c.i_max, c.ts * 1e6);
    printf("Design: wn=%.2f rad/s zeta=%.2f -> Kp=%.4f A/(rad/s)  Ki=%.3f A/rad\n",
           c.wn, c.zeta, g.kp, g.ki);
    printf("Q16.16: Kp=%d  Ki*Ts=%d (LSB)\n\n", q16_from_double(g.kp),
           q16_from_double(g.ki * c.ts));

    /* 1. small-signal step */
    int n = (int)(0.4 / c.ts);
    double *w = malloc(n * sizeof *w);
    loop_run(&c, &g, AW_CLAMP, W_BIAS, ref_small, NULL, n, w, NULL);
    printf("1) %.0f rad/s step at %.0f rad/s, sampled loop vs analytic continuous\n",
           DW_SMALL, W_BIAS);
    printf("   t [ms]   sim      analytic\n");
    double err = 0;
    for (int k = 0; k < n; k++) {
        double y = (w[k] - W_BIAS) / DW_SMALL;
        double ya = loop_analytic_step(&c, &g, k * c.ts);
        err = fmax(err, fabs(y - ya));
        if (k % 20 == 0 && k <= 160)
            printf("   %5d   %.4f   %.4f\n", k, y, ya);
    }
    printf("   max |sim - analytic| = %.4f of step\n\n", err);
    free(w);

    /* 2. large step, saturated */
    double w_inf = c.kt * c.i_max / c.b;
    double t_slew = c.J / c.b * log(w_inf / (w_inf - 0.9 * W_BIG));
    printf("2) 0 -> %.0f rad/s (current-limited, analytic t90 = %.3f s)\n",
           W_BIG, t_slew);
    big_step_metrics(&c, &g, AW_NONE, "no AW");
    big_step_metrics(&c, &g, AW_CLAMP, "clamp AW");
    printf("\n");

    /* 3. load torque step at constant speed */
    n = (int)(2.0 / c.ts);
    w = malloc(n * sizeof *w);
    loop_run(&c, &g, AW_CLAMP, W_BIG, ref_hold, load_step, n, w, NULL);
    double dip = 0;
    for (int k = 0; k < n; k++) dip = fmin(dip, w[k] - W_BIG);
    /* W/Td = -s / (J (s^2 + 2 zeta wn s + wn^2)); step response is
     * -(1/(J wd)) e^{-s t} sin(wd t), peaking at wd t = atan(wd/s). */
    double sg = c.zeta * c.wn, wd = c.wn * sqrt(1 - c.zeta * c.zeta);
    double tp = atan(wd / sg) / wd;
    double dip_a = -TAU_STEP / (c.J * wd) * exp(-sg * tp) * sin(wd * tp);
    printf("3) %.1f mN m load step at %.0f rad/s\n", TAU_STEP * 1e3, W_BIG);
    printf("   peak dip %.4f rad/s (analytic %.4f), final error %.2e rad/s\n",
           dip, dip_a, w[n - 1] - W_BIG);
    free(w);
    return 0;
}
