/*
 * Astron Built this AP
 * Project 017 | Attitude control | Rate- and torque-limited eigenaxis slew with reaction wheels
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-06 IST
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "scenario.h"

#define D2R (M_PI / 180.0)

static int fails = 0;

static void check(int ok, const char *what, double got, double lim)
{
    printf("[%s] %-58s %12.4e (limit %.3e)\n", ok ? "PASS" : "FAIL", what, got, lim);
    if (!ok) fails++;
}

/* --- angle history capture --- */
typedef struct { double *th; int n; long every; long i; } hist;

static void cap(double t, const sc_state *s, double ang, const vec3 *tau, void *p)
{
    (void)t; (void)s; (void)tau;
    hist *h = p;
    if (h->i % h->every == 0 && h->i / h->every < h->n)
        h->th[h->i / h->every] = ang;
    h->i++;
}

static const sim_opts OPTS = {.t_end = 600.0, .dt = 0.05, .settle_deg = 0.01, .settle_rate = 2e-5};

static void test_algebra(void)
{
    spacecraft sc = scenario_spacecraft();
    double err = 0.0;
    for (int j = 0; j < 3; j++) {
        vec3 ej = v3(j == 0, j == 1, j == 2);
        vec3 col = m3_mul(&sc.J, m3_mul(&sc.Jinv, ej));
        for (int i = 0; i < 3; i++)
            err = fmax(err, fabs(col.v[i] - (i == j)));
    }
    check(err < 1e-14, "J * Jinv = I", err, 1e-14);

    /* q_axis_angle then q_error_angle recovers the angle; rotate 90 deg about z. */
    quat a = q_axis_angle(v3(0.2, -0.4, 0.9), 0.3);
    quat b = q_mul(a, q_axis_angle(v3(1, 2, -1), 1.234));
    check(fabs(q_error_angle(a, b) - 1.234) < 1e-13, "relative angle of composed rotation",
          fabs(q_error_angle(a, b) - 1.234), 1e-13);
    vec3 r = q_rotate(q_axis_angle(v3(0, 0, 1), M_PI / 2), v3(1, 0, 0));
    double re = v3_norm(v3_add(r, v3(0, -1, 0)));
    check(re < 1e-15, "Rz(90) x_hat = y_hat", re, 1e-15);

    /* Bang-coast-bang: 1 rad, alpha 0.01, w_max 0.05 -> 1/0.05 + 0.05/0.01 = 25 s.
     * Triangular: 0.01 rad -> 2 sqrt(0.01/0.01) = 2 s. */
    check(fabs(slew_bcb_time(1.0, 0.01, 0.05) - 25.0) < 1e-12, "BCB time, trapezoidal", slew_bcb_time(1.0, 0.01, 0.05), 25.0);
    check(fabs(slew_bcb_time(0.01, 0.01, 0.05) - 2.0) < 1e-12, "BCB time, triangular", slew_bcb_time(0.01, 0.01, 0.05), 2.0);
}

static void test_design(void)
{
    spacecraft sc = scenario_spacecraft();
    vec3 h0 = scenario_h0();
    vec3 axis = v3(1, 2, -1);
    slew_gains g = slew_design(&sc, axis, v3_norm(h0), SCENARIO_ZETA, SCENARIO_WN,
                               SCENARIO_MARGIN, SCENARIO_BRAKE, 0.0);
    vec3 Ju = m3_mul(&sc.J, v3_unit(axis));
    double jmax = fmax(fabs(Ju.v[0]), fmax(fabs(Ju.v[1]), fabs(Ju.v[2])));
    /* Hand numbers: (J u)_y = (-15 + 3000 + 8)/sqrt(6) = 1221.9 kg m^2 */
    check(fabs(jmax - 2993.0 / sqrt(6.0)) < 1e-9, "dominant wheel lever |(J u)_y|", jmax, 2993.0 / sqrt(6.0));
    double wmax = (0.8 * 16.0 - v3_norm(h0)) / jmax;
    check(fabs(g.w_max - wmax) < 1e-15, "w_max = (0.8 h_max - |h0|) / |(J u)|max", g.w_max, wmax);
    check(fabs(g.alpha - 0.48 / jmax) < 1e-15, "alpha = 0.8 tau_max / |(J u)|max", g.alpha, 0.48 / jmax);
    /* Worst-case wheel torque: alpha*jmax + |w x H| <= tau_max with |H| = |h0|. */
    double worst = g.alpha * jmax + g.w_max * v3_norm(h0);
    check(worst <= sc.tau_max, "torque budget incl. gyroscopic term", worst, sc.tau_max);
    /* Rate command meets w_max where sqrt(4 a e) = w_max, and is continuous. */
    double e_star = g.w_max * g.w_max / (4 * g.a);
    check(fabs(slew_rate_command(&g, e_star) - g.w_max) < 1e-15, "rate command corner", slew_rate_command(&g, e_star), g.w_max);
}

/* Reference 90 deg slew: eigenaxis, limits, conservation, scalar match. */
static void test_reference_slew(void)
{
    spacecraft sc = scenario_spacecraft();
    vec3 h0 = scenario_h0();
    vec3 axis = v3(1, 2, -1);
    double ang = 90 * D2R;
    quat q0 = q_axis_angle(v3(0.3, -0.5, 0.8), 0.7);
    quat qc = q_mul(q0, q_axis_angle(axis, ang));
    slew_gains g = slew_design(&sc, axis, v3_norm(h0), SCENARIO_ZETA, SCENARIO_WN,
                               SCENARIO_MARGIN, SCENARIO_BRAKE, 0.0);
    sc_state s0 = {q0, v3(0, 0, 0), h0};

    enum { NS = 601 };
    static double th3[NS], th1[NS];
    hist h = {th3, NS, 20, 0};   /* 1 s samples */
    sim_result r = slew_simulate(CTRL_EIGENAXIS, &sc, &g, s0, qc, &OPTS, cap, &h);

    check(r.max_axis_dev < 1e-6, "eigenaxis: rotation-axis deviation (rad)", r.max_axis_dev, 1e-6);
    check(r.clip_time == 0.0 && r.max_tau <= 0.8 * sc.tau_max + v3_norm(h0) * g.w_max,
          "eigenaxis: wheel torque within design budget (N m)", r.max_tau, 0.8 * sc.tau_max + v3_norm(h0) * g.w_max);
    check(r.max_h <= 0.8 * sc.h_max, "eigenaxis: wheel momentum within 80% capacity (N m s)", r.max_h, 0.8 * sc.h_max);
    check(r.max_rate <= g.w_max * (1 + 1e-9), "eigenaxis: body rate never exceeds w_max (rad/s)", r.max_rate, g.w_max);
    check(r.max_rate >= 0.999 * g.w_max, "eigenaxis: coast phase reaches w_max (rad/s)", r.max_rate, 0.999 * g.w_max);
    check(r.H_drift < 1e-10, "inertial angular momentum conserved (relative)", r.H_drift, 1e-10);
    check(r.q_norm_err < 1e-12, "quaternion norm drift per step", r.q_norm_err, 1e-12);
    check(r.t_settle > 0 && r.final_angle < 0.01 * D2R, "settles inside 0.01 deg", r.final_angle, 0.01 * D2R);

    /* With the gyroscopic term cancelled, the 3D slew is exactly the scalar
     * eigenaxis equation. Integrate that independently with a 10x finer step. */
    slew_scalar_reference(&g, ang, 600.0, 0.005, 1.0, th1, NS);
    double dmax = 0.0;
    for (int i = 0; i < NS; i++)
        dmax = fmax(dmax, fabs(th3[i] - fabs(th1[i])));
    check(dmax < 1e-6, "3D angle vs scalar eigenaxis ODE, max |d theta| (rad)", dmax, 1e-6);

    /* Physical lower bound: with |theta''| <= alpha and |theta'| <= w_max no
     * profile reaches the 0.01 deg box sooner than bang-coast-bang. */
    double tb = slew_bcb_time(ang - 0.01 * D2R, g.alpha, g.w_max);
    check(r.t_settle >= tb, "settle time not below BCB bound (s)", r.t_settle, tb);
    check(r.t_settle <= 1.3 * slew_bcb_time(ang, g.alpha, g.w_max), "settle time within 30% of BCB (s)",
          r.t_settle, 1.3 * slew_bcb_time(ang, g.alpha, g.w_max));
}

/* Inside the linear region the slew is the textbook second-order response. */
static void test_small_angle(void)
{
    spacecraft sc = scenario_spacecraft();
    vec3 axis = v3(-0.4, 0.1, 0.9);
    slew_gains g = slew_design(&sc, axis, 0.0, SCENARIO_ZETA, SCENARIO_WN, SCENARIO_MARGIN, SCENARIO_BRAKE, 0.0);
    double th0 = 0.02;   /* k th0 = 2e-4 < alpha, |e| well inside a c^2/k^2 */
    check(g.k * th0 < g.alpha && 0.5 * th0 < g.a * g.c * g.c / (g.k * g.k),
          "small-angle case stays linear (k th0 / alpha)", g.k * th0 / g.alpha, 1.0);

    quat q0 = {0, 0, 0, 1};
    quat qc = q_axis_angle(axis, th0);
    sc_state s0 = {q0, v3(0, 0, 0), v3(0, 0, 0)};
    enum { NS = 121 };
    static double th[NS];
    hist h = {th, NS, 20, 0};
    sim_opts o = OPTS;
    o.t_end = 120.0;
    slew_simulate(CTRL_EIGENAXIS, &sc, &g, s0, qc, &o, cap, &h);

    double z = g.zeta, wn = g.wn, wd = wn * sqrt(1 - z * z), err = 0.0;
    for (int i = 0; i < NS; i++) {
        double t = i;
        double an = th0 * exp(-z * wn * t) * (cos(wd * t) + z / sqrt(1 - z * z) * sin(wd * t));
        err = fmax(err, fabs(th[i] - fabs(an)));
    }
    check(err < 1e-4 * th0, "small slew vs analytic 2nd-order response (rad)", err, 1e-4 * th0);
}

/* Without the gyroscopic term the rotation leaves the eigenaxis, and the
 * unlimited scalar PD saturates wheels: the reasons for the full law. */
static void test_baselines(void)
{
    spacecraft sc = scenario_spacecraft();
    vec3 h0 = scenario_h0();
    vec3 axis = v3(1, 2, -1);
    quat q0 = q_axis_angle(v3(0.3, -0.5, 0.8), 0.7);
    quat qc = q_mul(q0, q_axis_angle(axis, 90 * D2R));
    slew_gains g = slew_design(&sc, axis, v3_norm(h0), SCENARIO_ZETA, SCENARIO_WN,
                               SCENARIO_MARGIN, SCENARIO_BRAKE, 0.0);
    sc_state s0 = {q0, v3(0, 0, 0), h0};
    sim_result a = slew_simulate(CTRL_EIGENAXIS, &sc, &g, s0, qc, &OPTS, NULL, NULL);
    sim_result b = slew_simulate(CTRL_EIGENAXIS_NOGYRO, &sc, &g, s0, qc, &OPTS, NULL, NULL);
    sim_result c = slew_simulate(CTRL_NAIVE_PD, &sc, &g, s0, qc, &OPTS, NULL, NULL);
    check(b.max_axis_dev > 1e4 * a.max_axis_dev, "no gyro comp: axis deviation >> compensated (rad)", b.max_axis_dev, 1e4 * a.max_axis_dev);
    check(c.clip_time > 10.0, "scalar PD: wheel torque saturates (s clipped)", c.clip_time, 10.0);
    check(c.max_h > sc.h_max, "scalar PD: exceeds wheel momentum capacity (N m s)", c.max_h, sc.h_max);
}

/* Random axes, angles, bias momenta and inertia-frame attitudes. */
static void test_monte_carlo(void)
{
    spacecraft sc = scenario_spacecraft();
    srand(17);
    int ok_axis = 1, ok_lim = 1, ok_settle = 1, ok_bound = 1;
    double worst_dev = 0, worst_ratio = 0;
    for (int n = 0; n < 40; n++) {
        vec3 axis, h0, qa;
        for (int i = 0; i < 3; i++) {
            axis.v[i] = 2.0 * rand() / RAND_MAX - 1.0;
            h0.v[i] = 3.0 * rand() / RAND_MAX - 1.5;
            qa.v[i] = 2.0 * rand() / RAND_MAX - 1.0;
        }
        double ang = (2.0 + 176.0 * rand() / RAND_MAX) * D2R;
        quat q0 = q_axis_angle(qa, 3.0 * rand() / RAND_MAX);
        quat qc = q_mul(q0, q_axis_angle(axis, ang));
        slew_gains g = slew_design(&sc, axis, v3_norm(h0), SCENARIO_ZETA, SCENARIO_WN,
                                   SCENARIO_MARGIN, SCENARIO_BRAKE, 0.0);
        sc_state s0 = {q0, v3(0, 0, 0), h0};
        sim_opts o = OPTS;
        o.t_end = 1500.0;
        sim_result r = slew_simulate(CTRL_EIGENAXIS, &sc, &g, s0, qc, &o, NULL, NULL);
        worst_dev = fmax(worst_dev, r.max_axis_dev);
        if (r.max_axis_dev > 1e-6) ok_axis = 0;
        if (r.clip_time > 0 || r.max_h > 0.8 * sc.h_max || r.max_rate > g.w_max * (1 + 1e-9)) ok_lim = 0;
        if (r.t_settle < 0) ok_settle = 0;
        double tb = slew_bcb_time(ang - 0.01 * D2R, g.alpha, g.w_max);
        if (r.t_settle < tb) ok_bound = 0;
        worst_ratio = fmax(worst_ratio, r.t_settle / slew_bcb_time(ang, g.alpha, g.w_max));
    }
    check(ok_axis, "MC 40 slews: eigenaxis held (worst rad)", worst_dev, 1e-6);
    check(ok_lim, "MC 40 slews: no torque clip, |h| <= 0.8 h_max, |w| <= w_max", 0, 0);
    check(ok_settle, "MC 40 slews: all settle", 0, 0);
    check(ok_bound, "MC 40 slews: none beats BCB bound (worst T/T_bcb)", worst_ratio, 0);
}

int main(void)
{
    test_algebra();
    test_design();
    test_reference_slew();
    test_small_angle();
    test_baselines();
    test_monte_carlo();
    printf("\n%s (%d failure%s)\n", fails ? "FAILED" : "ALL PASSED", fails, fails == 1 ? "" : "s");
    return fails ? 1 : 0;
}
