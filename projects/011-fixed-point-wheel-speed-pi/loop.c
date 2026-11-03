/*
 * Astron Built this AP
 * Project 011 | Embedded firmware | Fixed-point PI speed loop for a reaction wheel with anti-windup
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-03 IST
 */
#include "loop.h"

#include <math.h>
#include <stddef.h>

/* A 0.2 Nms class small-satellite wheel: 2e-4 kg m^2 rotor, 20 mNm/A motor,
 * 1.5 A driver limit -> 30 mNm max torque. 1 kHz control loop. */
loop_cfg_t loop_default_cfg(void)
{
    loop_cfg_t c = {
        .J = 2.0e-4, .kt = 0.02, .b = 1.0e-6,
        .i_max = 1.5, .ts = 1.0e-3,
        .wn = 2.0 * M_PI * 8.0, .zeta = 0.8,
    };
    return c;
}

/* Closed-loop characteristic polynomial with u = Kp e + Ki int(e):
 *   J s^2 + (b + Kt Kp) s + Kt Ki = J (s^2 + 2 zeta wn s + wn^2)
 * so Kp = (2 zeta wn J - b) / Kt and Ki = wn^2 J / Kt. */
pi_gains_t loop_design(const loop_cfg_t *c)
{
    pi_gains_t g;
    g.kp = (2.0 * c->zeta * c->wn * c->J - c->b) / c->kt;
    g.ki = c->wn * c->wn * c->J / c->kt;
    return g;
}

/* W/R = (a s + wn^2) / (s^2 + 2 zeta wn s + wn^2),  a = Kt Kp / J.
 * Step response = standard 2nd-order step + a * (its impulse response):
 *   y = 1 - e^{-st}(cos wd t + s/wd sin wd t) + (a/wd) e^{-st} sin wd t
 * with s = zeta wn, wd = wn sqrt(1 - zeta^2). Valid for zeta < 1. */
double loop_analytic_step(const loop_cfg_t *c, const pi_gains_t *g, double t)
{
    double a = c->kt * g->kp / c->J;
    double sg = c->zeta * c->wn;
    double wd = c->wn * sqrt(1.0 - c->zeta * c->zeta);
    double ex = exp(-sg * t);
    return 1.0 - ex * (cos(wd * t) + sg / wd * sin(wd * t))
           + a / wd * ex * sin(wd * t);
}

void loop_run(const loop_cfg_t *cfg, const pi_gains_t *g, aw_mode_t aw,
              double w0, profile_fn ref, profile_fn load, int n,
              double *w_out, double *i_out)
{
    pi_q16_t pi;
    pi_q16_init(&pi, g->kp, g->ki, cfg->ts, cfg->i_max, aw);
    wheel_t wh = { .J = cfg->J, .kt = cfg->kt, .b = cfg->b, .w = w0 };

    /* Start from the integrator value that holds w0 against drag, so a
     * nonzero initial speed is an equilibrium rather than a transient. */
    pi.integ = q16_from_double(cfg->b * w0 / cfg->kt);

    for (int k = 0; k < n; k++) {
        double t = k * cfg->ts;
        if (w_out) w_out[k] = wh.w;
        q16_t u = pi_q16_step(&pi, q16_from_double(ref(t)),
                              q16_from_double(wh.w));
        double i = q16_to_double(u);
        if (i_out) i_out[k] = i;
        wheel_step(&wh, i, load ? load(t) : 0.0, cfg->ts);
    }
}

void pi_dbl_init(pi_dbl_t *c, double kp, double ki, double ts, double u_max,
                 aw_mode_t aw)
{
    c->kp = kp; c->ki_dt = ki * ts; c->u_max = u_max; c->integ = 0; c->aw = aw;
}

double pi_dbl_step(pi_dbl_t *c, double ref, double meas)
{
    double e = ref - meas;
    double u_raw = c->kp * e + c->integ;
    double u = fmin(fmax(u_raw, -c->u_max), c->u_max);
    double di = c->ki_dt * e;
    if (c->aw == AW_CLAMP) {
        int hi = u_raw >= c->u_max && e > 0;
        int lo = u_raw <= -c->u_max && e < 0;
        if (!hi && !lo)
            c->integ = fmin(fmax(c->integ + di, -c->u_max), c->u_max);
    } else {
        c->integ += di;
    }
    return u;
}
