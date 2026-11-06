/*
 * Astron Built this AP
 * Project 017 | Attitude control | Rate- and torque-limited eigenaxis slew with reaction wheels
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-06 IST
 */
#include "slew.h"

#include <math.h>
#include <stddef.h>

/* ---------- small linear algebra ---------- */

vec3 v3(double x, double y, double z) { vec3 r = {{x, y, z}}; return r; }

double v3_dot(vec3 a, vec3 b)
{
    return a.v[0] * b.v[0] + a.v[1] * b.v[1] + a.v[2] * b.v[2];
}

vec3 v3_cross(vec3 a, vec3 b)
{
    return v3(a.v[1] * b.v[2] - a.v[2] * b.v[1],
              a.v[2] * b.v[0] - a.v[0] * b.v[2],
              a.v[0] * b.v[1] - a.v[1] * b.v[0]);
}

double v3_norm(vec3 a) { return sqrt(v3_dot(a, a)); }

vec3 v3_scale(vec3 a, double s) { return v3(a.v[0] * s, a.v[1] * s, a.v[2] * s); }

vec3 v3_add(vec3 a, vec3 b) { return v3(a.v[0] + b.v[0], a.v[1] + b.v[1], a.v[2] + b.v[2]); }

vec3 v3_unit(vec3 a) { return v3_scale(a, 1.0 / v3_norm(a)); }

vec3 m3_mul(const mat3 *A, vec3 x)
{
    vec3 r;
    for (int i = 0; i < 3; i++)
        r.v[i] = A->m[i][0] * x.v[0] + A->m[i][1] * x.v[1] + A->m[i][2] * x.v[2];
    return r;
}

bool m3_inverse(const mat3 *A, mat3 *out)
{
    const double (*a)[3] = A->m;
    double c00 = a[1][1] * a[2][2] - a[1][2] * a[2][1];
    double c01 = a[1][2] * a[2][0] - a[1][0] * a[2][2];
    double c02 = a[1][0] * a[2][1] - a[1][1] * a[2][0];
    double det = a[0][0] * c00 + a[0][1] * c01 + a[0][2] * c02;
    if (fabs(det) < 1e-300)
        return false;
    double id = 1.0 / det;
    out->m[0][0] = c00 * id;
    out->m[1][0] = c01 * id;
    out->m[2][0] = c02 * id;
    out->m[0][1] = (a[0][2] * a[2][1] - a[0][1] * a[2][2]) * id;
    out->m[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) * id;
    out->m[2][1] = (a[0][1] * a[2][0] - a[0][0] * a[2][1]) * id;
    out->m[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) * id;
    out->m[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) * id;
    out->m[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) * id;
    return true;
}

/* Hamilton product, scalar last. */
quat q_mul(quat a, quat b)
{
    quat r;
    r.w = a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z;
    r.x = a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y;
    r.y = a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x;
    r.z = a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w;
    return r;
}

quat q_conj(quat a) { quat r = {-a.x, -a.y, -a.z, a.w}; return r; }

quat q_axis_angle(vec3 axis, double angle)
{
    vec3 u = v3_unit(axis);
    double s = sin(0.5 * angle);
    quat r = {u.v[0] * s, u.v[1] * s, u.v[2] * s, cos(0.5 * angle)};
    return r;
}

vec3 q_rotate(quat q, vec3 v)
{
    quat p = {v.v[0], v.v[1], v.v[2], 0.0};
    quat r = q_mul(q_mul(q, p), q_conj(q));
    return v3(r.x, r.y, r.z);
}

static quat q_normalize(quat q)
{
    double n = sqrt(q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w);
    quat r = {q.x / n, q.y / n, q.z / n, q.w / n};
    return r;
}

/* Error quaternion q_cmd^* (x) q. Its vector part is the rotation axis times
 * sin(theta/2), identical in body and command frames. Sign is chosen so the
 * scalar part is non-negative, which selects the short way round. */
static quat q_error(quat q, quat q_cmd)
{
    quat e = q_mul(q_conj(q_cmd), q);
    if (e.w < 0.0) {
        e.x = -e.x; e.y = -e.y; e.z = -e.z; e.w = -e.w;
    }
    return e;
}

double q_error_angle(quat q, quat q_cmd)
{
    quat e = q_error(q, q_cmd);
    return 2.0 * atan2(sqrt(e.x * e.x + e.y * e.y + e.z * e.z), e.w);
}

/* ---------- design ---------- */

slew_gains slew_design(const spacecraft *sc, vec3 axis_body, double h0_norm,
                       double zeta, double wn, double margin, double brake_frac,
                       double rate_cap)
{
    /* For a rotation about unit axis u at rate w, wheel i must absorb
     * w*|(J u)_i| plus whatever bias momentum it already holds (at most |h0|
     * in any body direction). The same column bounds torque: accelerating at
     * alpha needs alpha*|(J u)_i| from wheel i. */
    vec3 u = v3_unit(axis_body);
    vec3 Ju = m3_mul(&sc->J, u);
    double jmax = 0.0;
    for (int i = 0; i < 3; i++)
        if (fabs(Ju.v[i]) > jmax)
            jmax = fabs(Ju.v[i]);

    slew_gains g;
    g.w_max = (margin * sc->h_max - h0_norm) / jmax;
    if (rate_cap > 0.0 && rate_cap < g.w_max)
        g.w_max = rate_cap;
    g.alpha = margin * sc->tau_max / jmax;

    /* Small angles: e ~ theta/2, so theta'' + c theta' + k theta = 0. */
    g.zeta = zeta;
    g.wn = wn;
    g.k = wn * wn;
    g.c = 2.0 * zeta * wn;
    g.a = brake_frac * g.alpha;
    return g;
}

double slew_rate_command(const slew_gains *g, double en)
{
    /* |e| = sin(theta/2) ~ theta/2, so sqrt(4 a |e|) ~ sqrt(2 a theta): the
     * rate from which a constant deceleration a stops exactly on target. */
    double w = sqrt(4.0 * g->a * en);
    return w < g->w_max ? w : g->w_max;
}

/* Saturation limit L(e) on |e|; inside it the law is the linear PD. */
static double sat_limit(const slew_gains *g, double en)
{
    return g->c / (2.0 * g->k) * slew_rate_command(g, en);
}

/* ---------- control laws ---------- */

vec3 slew_torque(ctrl_mode mode, const spacecraft *sc, const slew_gains *g,
                 const sc_state *s, quat q_cmd)
{
    quat qe = q_error(s->q, q_cmd);
    vec3 e = v3(qe.x, qe.y, qe.z);

    if (mode == CTRL_NAIVE_PD) {
        /* Same small-angle bandwidth, but a scalar inertia and no limits:
         * the textbook law that works for small corrections. */
        double jbar = (sc->J.m[0][0] + sc->J.m[1][1] + sc->J.m[2][2]) / 3.0;
        return v3_add(v3_scale(e, -2.0 * g->k * jbar), v3_scale(s->w, -g->c * jbar));
    }

    /* Saturating whole vectors, not components, keeps the commanded
     * acceleration parallel to e; with the gyroscopic term cancelled, w and e
     * then stay on the initial eigenaxis for the whole manoeuvre.
     * When sat_L is active the law reads w_dot = -c (w - w_cmd) with
     * w_cmd = -e_hat * min(sqrt(4a|e|), w_max): a rate loop riding the
     * braking curve. */
    double en = v3_norm(e);
    double L = sat_limit(g, en);
    vec3 es = (en > L) ? v3_scale(e, L / en) : e;
    vec3 acc = v3_add(v3_scale(es, -2.0 * g->k), v3_scale(s->w, -g->c));
    double an = v3_norm(acc);
    if (an > g->alpha)
        acc = v3_scale(acc, g->alpha / an);
    vec3 tau = m3_mul(&sc->J, acc);
    if (mode == CTRL_EIGENAXIS) {
        vec3 H = v3_add(m3_mul(&sc->J, s->w), s->h);
        tau = v3_add(tau, v3_cross(s->w, H));
    }
    return tau;
}

/* ---------- dynamics ---------- */

typedef struct { quat dq; vec3 dw; vec3 dh; } sc_deriv;

typedef struct {
    ctrl_mode mode;
    const spacecraft *sc;
    const slew_gains *g;
    quat q_cmd;
} sim_ctx;

static vec3 clip_torque(vec3 tau, double tmax, bool *clipped)
{
    *clipped = false;
    for (int i = 0; i < 3; i++) {
        if (tau.v[i] > tmax) { tau.v[i] = tmax; *clipped = true; }
        if (tau.v[i] < -tmax) { tau.v[i] = -tmax; *clipped = true; }
    }
    return tau;
}

/* J w_dot = tau - w x (J w + h),  h_dot = -tau,  q_dot = 1/2 q (x) [w, 0].
 * Control is evaluated continuously inside the derivative. */
static sc_deriv deriv(const sim_ctx *c, const sc_state *s)
{
    bool clipped;
    vec3 tau = clip_torque(slew_torque(c->mode, c->sc, c->g, s, c->q_cmd),
                           c->sc->tau_max, &clipped);
    vec3 H = v3_add(m3_mul(&c->sc->J, s->w), s->h);
    sc_deriv d;
    d.dw = m3_mul(&c->sc->Jinv, v3_add(tau, v3_scale(v3_cross(s->w, H), -1.0)));
    d.dh = v3_scale(tau, -1.0);
    quat wq = {s->w.v[0], s->w.v[1], s->w.v[2], 0.0};
    quat qd = q_mul(s->q, wq);
    d.dq.x = 0.5 * qd.x; d.dq.y = 0.5 * qd.y; d.dq.z = 0.5 * qd.z; d.dq.w = 0.5 * qd.w;
    return d;
}

static sc_state step_add(const sc_state *s, const sc_deriv *d, double h)
{
    sc_state r;
    r.q.x = s->q.x + h * d->dq.x;
    r.q.y = s->q.y + h * d->dq.y;
    r.q.z = s->q.z + h * d->dq.z;
    r.q.w = s->q.w + h * d->dq.w;
    r.w = v3_add(s->w, v3_scale(d->dw, h));
    r.h = v3_add(s->h, v3_scale(d->dh, h));
    return r;
}

static sc_state rk4(const sim_ctx *c, const sc_state *s, double h)
{
    sc_deriv k1 = deriv(c, s);
    sc_state s2 = step_add(s, &k1, 0.5 * h);
    sc_deriv k2 = deriv(c, &s2);
    sc_state s3 = step_add(s, &k2, 0.5 * h);
    sc_deriv k3 = deriv(c, &s3);
    sc_state s4 = step_add(s, &k3, h);
    sc_deriv k4 = deriv(c, &s4);

    sc_deriv k;
    k.dq.x = (k1.dq.x + 2 * k2.dq.x + 2 * k3.dq.x + k4.dq.x) / 6.0;
    k.dq.y = (k1.dq.y + 2 * k2.dq.y + 2 * k3.dq.y + k4.dq.y) / 6.0;
    k.dq.z = (k1.dq.z + 2 * k2.dq.z + 2 * k3.dq.z + k4.dq.z) / 6.0;
    k.dq.w = (k1.dq.w + 2 * k2.dq.w + 2 * k3.dq.w + k4.dq.w) / 6.0;
    for (int i = 0; i < 3; i++) {
        k.dw.v[i] = (k1.dw.v[i] + 2 * k2.dw.v[i] + 2 * k3.dw.v[i] + k4.dw.v[i]) / 6.0;
        k.dh.v[i] = (k1.dh.v[i] + 2 * k2.dh.v[i] + 2 * k3.dh.v[i] + k4.dh.v[i]) / 6.0;
    }
    return step_add(s, &k, h);
}

sim_result slew_simulate(ctrl_mode mode, const spacecraft *sc, const slew_gains *g,
                         sc_state s0, quat q_cmd, const sim_opts *o,
                         sim_logger log, void *ctx)
{
    sim_ctx c = {mode, sc, g, q_cmd};
    sim_result r = {0};
    sc_state s = s0;

    quat qe0 = q_error(s0.q, q_cmd);
    vec3 axis0 = v3_unit(v3(qe0.x, qe0.y, qe0.z));
    vec3 H0 = q_rotate(s0.q, v3_add(m3_mul(&sc->J, s0.w), s0.h));
    double H0n = v3_norm(H0);
    double settle_rad = o->settle_deg * M_PI / 180.0;
    double last_outside = 0.0;

    long n = lround(o->t_end / o->dt);
    for (long i = 0; i <= n; i++) {
        double t = i * o->dt;
        bool clipped;
        vec3 tau_raw = slew_torque(mode, sc, g, &s, q_cmd);
        vec3 tau = clip_torque(tau_raw, sc->tau_max, &clipped);
        double ang = q_error_angle(s.q, q_cmd);

        for (int j = 0; j < 3; j++) {
            if (fabs(tau_raw.v[j]) > r.max_tau) r.max_tau = fabs(tau_raw.v[j]);
            if (fabs(s.h.v[j]) > r.max_h) r.max_h = fabs(s.h.v[j]);
        }
        double wn = v3_norm(s.w);
        if (wn > r.max_rate) r.max_rate = wn;
        if (clipped && i < n) r.clip_time += o->dt;

        /* Rotation-axis deviation, only while the axis is well defined. */
        quat qe = q_error(s.q, q_cmd);
        vec3 e = v3(qe.x, qe.y, qe.z);
        if (v3_norm(e) > 0.01 * v3_norm(v3(qe0.x, qe0.y, qe0.z))) {
            double cs = fabs(v3_dot(v3_unit(e), axis0));
            double dev = acos(cs > 1.0 ? 1.0 : cs);
            if (dev > r.max_axis_dev) r.max_axis_dev = dev;
        }

        vec3 H = q_rotate(s.q, v3_add(m3_mul(&sc->J, s.w), s.h));
        double dH = v3_norm(v3_add(H, v3_scale(H0, -1.0))) / (H0n > 1e-12 ? H0n : 1.0);
        if (dH > r.H_drift) r.H_drift = dH;

        if (ang > settle_rad || wn > o->settle_rate)
            last_outside = t;

        if (log)
            log(t, &s, ang, &tau, ctx);
        if (i == n) {
            r.final_angle = ang;
            break;
        }
        s = rk4(&c, &s, o->dt);
        /* RK4 does not preserve |q| exactly; record the per-step drift as an
         * integration-accuracy check, then renormalise. */
        double qn = sqrt(s.q.x * s.q.x + s.q.y * s.q.y + s.q.z * s.q.z + s.q.w * s.q.w);
        if (fabs(qn - 1.0) > r.q_norm_err) r.q_norm_err = fabs(qn - 1.0);
        s.q = q_normalize(s.q);
    }
    r.t_settle = (last_outside < o->t_end) ? last_outside : -1.0;
    return r;
}

/* ---------- references ---------- */

static void scalar_f(const slew_gains *g, double th, double thd, double *a)
{
    double e = sin(0.5 * th);
    double L = sat_limit(g, fabs(e));
    if (fabs(e) > L) e = copysign(L, e);
    double acc = -2.0 * g->k * e - g->c * thd;
    if (fabs(acc) > g->alpha) acc = copysign(g->alpha, acc);
    *a = acc;
}

void slew_scalar_reference(const slew_gains *g, double th0, double t_end,
                           double dt, double out_dt, double *theta, int n_out)
{
    double th = th0, thd = 0.0;
    long n = lround(t_end / dt);
    long every = lround(out_dt / dt);
    for (long i = 0; i <= n; i++) {
        if (i % every == 0) {
            long j = i / every;
            if (j < n_out) theta[j] = th;
        }
        double a1, a2, a3, a4;
        scalar_f(g, th, thd, &a1);
        scalar_f(g, th + 0.5 * dt * thd, thd + 0.5 * dt * a1, &a2);
        scalar_f(g, th + 0.5 * dt * (thd + 0.5 * dt * a1), thd + 0.5 * dt * a2, &a3);
        scalar_f(g, th + dt * (thd + 0.5 * dt * a2), thd + dt * a3, &a4);
        th += dt * (thd + dt * (a1 + a2 + a3) / 6.0);
        thd += dt * (a1 + 2 * a2 + 2 * a3 + a4) / 6.0;
    }
}

double slew_bcb_time(double angle, double alpha, double w_max)
{
    /* Triangular profile if the rate limit is never reached. */
    if (angle <= w_max * w_max / alpha)
        return 2.0 * sqrt(angle / alpha);
    return angle / w_max + w_max / alpha;
}
