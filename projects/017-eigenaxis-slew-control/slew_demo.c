/*
 * Astron Built this AP
 * Project 017 | Attitude control | Rate- and torque-limited eigenaxis slew with reaction wheels
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-06 IST
 *
 * 90 deg retarget of an agile imager. Compares the saturated eigenaxis law
 * against the same law without gyroscopic compensation and against a plain
 * scalar-gain PD, then sweeps slew angle against the bang-coast-bang bound.
 * Usage: ./slew_demo [history.csv]
 */
#include <math.h>
#include <stdio.h>

#include "scenario.h"

#define D2R (M_PI / 180.0)
#define R2D (180.0 / M_PI)

typedef struct { FILE *f; long n; } csv_ctx;

static void csv_log(double t, const sc_state *s, double ang, const vec3 *tau, void *p)
{
    csv_ctx *c = p;
    if (c->n++ % 20)   /* 1 s rows at dt = 0.05 s */
        return;
    fprintf(c->f, "%.2f,%.6f,%.6e,%.6e,%.6e,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n",
            t, ang * R2D, s->w.v[0], s->w.v[1], s->w.v[2],
            tau->v[0], tau->v[1], tau->v[2], s->h.v[0], s->h.v[1], s->h.v[2]);
}

int main(int argc, char **argv)
{
    spacecraft sc = scenario_spacecraft();
    vec3 h0 = scenario_h0();
    vec3 axis = v3(1.0, 2.0, -1.0);
    double angle = 90.0 * D2R;

    quat q0 = q_axis_angle(v3(0.3, -0.5, 0.8), 0.7);
    quat q_cmd = q_mul(q0, q_axis_angle(axis, angle));

    slew_gains g = slew_design(&sc, axis, v3_norm(h0), SCENARIO_ZETA, SCENARIO_WN,
                               SCENARIO_MARGIN, SCENARIO_BRAKE, 0.0);
    printf("Agile imager, 3 body-axis wheels: tau_max = %.2f N m, h_max = %.1f N m s\n",
           sc.tau_max, sc.h_max);
    printf("Bias wheel momentum |h0| = %.2f N m s, slew %.1f deg about body [1 2 -1]\n\n",
           v3_norm(h0), angle * R2D);
    printf("Design (margin %.0f%%)\n", SCENARIO_MARGIN * 100);
    printf("  w_max = %.4f deg/s   alpha = %.3e rad/s^2   brake a = %.3e rad/s^2\n",
           g.w_max * R2D, g.alpha, g.a);
    printf("  k = %.4e  c = %.4e  (wn = %.3f rad/s, zeta = %.2f)\n", g.k, g.c, g.wn, g.zeta);
    /* Linear region: |e| < L(e)  <=>  |e| < a c^2 / k^2. */
    printf("  linear PD region below %.2f deg error\n\n",
           2.0 * asin(g.a * g.c * g.c / (g.k * g.k)) * R2D);

    sim_opts o = {.t_end = 600.0, .dt = 0.05, .settle_deg = 0.01, .settle_rate = 2e-5};
    sc_state s0 = {q0, v3(0, 0, 0), h0};

    const char *names[] = {"eigenaxis + gyro comp", "eigenaxis, no gyro comp", "scalar PD, unlimited"};
    ctrl_mode modes[] = {CTRL_EIGENAXIS, CTRL_EIGENAXIS_NOGYRO, CTRL_NAIVE_PD};
    printf("%-24s %9s %10s %9s %9s %9s %9s\n", "controller", "settle s", "axis dev", "tau pk",
           "clip s", "|h| pk", "rate pk");
    printf("%-24s %9s %10s %9s %9s %9s %9s\n", "", "", "deg", "N m", "", "N m s", "deg/s");
    for (int m = 0; m < 3; m++) {
        csv_ctx c = {0};
        if (m == 0 && argc > 1) {
            c.f = fopen(argv[1], "w");
            if (c.f)
                fprintf(c.f, "t,angle_deg,wx,wy,wz,tx,ty,tz,hx,hy,hz\n");
        }
        sim_result r = slew_simulate(modes[m], &sc, &g, s0, q_cmd, &o,
                                     c.f ? csv_log : NULL, &c);
        if (c.f)
            fclose(c.f);
        char ts[16];
        if (r.t_settle < 0) snprintf(ts, sizeof ts, "never");
        else snprintf(ts, sizeof ts, "%.1f", r.t_settle);
        printf("%-24s %9s %10.2e %9.3f %9.1f %9.2f %9.4f\n", names[m], ts,
               r.max_axis_dev * R2D, r.max_tau, r.clip_time, r.max_h, r.max_rate * R2D);
    }
    printf("\nBang-coast-bang bound with the same alpha and w_max: %.1f s\n\n",
           slew_bcb_time(angle, g.alpha, g.w_max));

    printf("Angle sweep about the same axis (eigenaxis law)\n");
    printf("%8s %10s %10s %7s\n", "deg", "settle s", "BCB s", "ratio");
    double sweep[] = {5, 15, 30, 60, 90, 120, 150, 175};
    for (int i = 0; i < 8; i++) {
        double a = sweep[i] * D2R;
        quat qc = q_mul(q0, q_axis_angle(axis, a));
        sim_opts os = o;
        os.t_end = 900.0;
        sim_result r = slew_simulate(CTRL_EIGENAXIS, &sc, &g, s0, qc, &os, NULL, NULL);
        double tb = slew_bcb_time(a, g.alpha, g.w_max);
        printf("%8.0f %10.1f %10.1f %7.2f\n", sweep[i], r.t_settle, tb, r.t_settle / tb);
    }
    return 0;
}
