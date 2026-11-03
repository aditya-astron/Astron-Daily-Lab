/*
 * Astron Built this AP
 * Project 011 | Embedded firmware | Fixed-point PI speed loop for a reaction wheel with anti-windup
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-03 IST
 */
#include "pi_q16.h"

#include <math.h>

q16_t q16_from_double(double x) { return (q16_t)lround(x * Q16_ONE); }

double q16_to_double(q16_t x) { return (double)x / Q16_ONE; }

/* Product in int64 then round-half-up back to Q16.16. Truncation here would
 * bias the integrator by up to 1 LSB per sample, which over thousands of
 * samples shows up as a slow steady-state drift. */
q16_t q16_mul(q16_t a, q16_t b)
{
    int64_t p = (int64_t)a * (int64_t)b;
    return (q16_t)((p + (1 << 15)) >> 16);
}

static q16_t sat(q16_t x, q16_t lim)
{
    if (x > lim) return lim;
    if (x < -lim) return -lim;
    return x;
}

void pi_q16_init(pi_q16_t *c, double kp, double ki, double ts,
                 double u_max, aw_mode_t aw)
{
    c->kp = q16_from_double(kp);
    c->ki_dt = q16_from_double(ki * ts);
    c->u_max = q16_from_double(u_max);
    c->integ = 0;
    c->aw = aw;
}

/* u = Kp*e + I,  I[k+1] = I[k] + Ki*Ts*e  (backward-rectangle integrator).
 *
 * Anti-windup by conditional integration: when the output is pinned at the
 * limit and the error would push it further in, the integrator is frozen.
 * The integrator is also bounded to +/-u_max on its own, since a larger
 * value can never be needed to hold any reachable speed. */
q16_t pi_q16_step(pi_q16_t *c, q16_t ref, q16_t meas)
{
    q16_t e = ref - meas;
    q16_t p = q16_mul(c->kp, e);
    int64_t u_raw = (int64_t)p + c->integ;
    q16_t u = (q16_t)(u_raw > c->u_max ? c->u_max
                      : u_raw < -c->u_max ? -c->u_max : u_raw);

    q16_t di = q16_mul(c->ki_dt, e);
    if (c->aw == AW_CLAMP) {
        int saturated_hi = u_raw >= c->u_max && e > 0;
        int saturated_lo = u_raw <= -c->u_max && e < 0;
        if (!saturated_hi && !saturated_lo)
            c->integ = sat(c->integ + di, c->u_max);
    } else {
        /* Unprotected integrator; int64 add then clip only at the Q16.16
         * range so the windup is physically what a naive loop would do. */
        int64_t s = (int64_t)c->integ + di;
        if (s > INT32_MAX) s = INT32_MAX;
        if (s < INT32_MIN) s = INT32_MIN;
        c->integ = (q16_t)s;
    }
    return u;
}
