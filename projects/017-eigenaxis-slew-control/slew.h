/*
 * Astron Built this AP
 * Project 017 | Attitude control | Rate- and torque-limited eigenaxis slew with reaction wheels
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-06 IST
 */
#ifndef SLEW_H
#define SLEW_H

#include <stdbool.h>

typedef struct { double v[3]; } vec3;
typedef struct { double m[3][3]; } mat3;
/* Scalar-last quaternion (x, y, z, w), body -> inertial. */
typedef struct { double x, y, z, w; } quat;

/* Spacecraft with three orthogonal reaction wheels on the body axes. */
typedef struct {
    mat3 J;            /* body inertia including wheels, kg m^2 */
    mat3 Jinv;
    double tau_max;    /* per-wheel torque limit, N m */
    double h_max;      /* per-wheel momentum capacity, N m s */
} spacecraft;

/* Gains of the commanded body acceleration (after gyroscopic cancellation)
 *   w_dot = sat_alpha( -2 k sat_L(e) - c w ),
 *   L(e)  = (c / 2k) min( sqrt(4 a |e|), w_max ),
 * where both saturations scale the whole vector, so direction is preserved. */
typedef struct {
    double k, c;
    double w_max;      /* slew rate limit, rad/s */
    double alpha;      /* angular acceleration limit, rad/s^2 */
    double a;          /* braking deceleration that shapes the rate command */
    double zeta, wn;   /* small-angle closed-loop damping and natural freq */
} slew_gains;

typedef enum {
    CTRL_EIGENAXIS,    /* J-scaled, saturated, gyroscopic compensation */
    CTRL_EIGENAXIS_NOGYRO,
    CTRL_NAIVE_PD      /* scalar-gain PD, no saturation, no compensation */
} ctrl_mode;

typedef struct {
    quat q;            /* attitude */
    vec3 w;            /* body rate, rad/s */
    vec3 h;            /* wheel momentum in body axes, N m s */
} sc_state;

typedef struct {
    double t_end, dt;
    double settle_deg;     /* settled when angle error stays below this */
    double settle_rate;    /* and |w| below this, rad/s */
} sim_opts;

typedef struct {
    double t_settle;       /* first time after which the settle box holds; <0 if never */
    double max_axis_dev;   /* max angle between w and the initial eigenaxis, rad */
    double max_tau;        /* max per-wheel commanded torque before clipping */
    double max_h;          /* max per-wheel momentum magnitude */
    double max_rate;       /* max |w| */
    double clip_time;      /* seconds with at least one wheel torque-saturated */
    double H_drift;        /* max relative drift of inertial angular momentum */
    double q_norm_err;     /* max | |q| - 1 | */
    double final_angle;    /* error angle at t_end, rad */
} sim_result;

/* Optional per-step logging callback. */
typedef void (*sim_logger)(double t, const sc_state *s, double angle, const vec3 *tau, void *ctx);

/* math helpers */
vec3 v3(double x, double y, double z);
double v3_dot(vec3 a, vec3 b);
vec3 v3_cross(vec3 a, vec3 b);
double v3_norm(vec3 a);
vec3 v3_scale(vec3 a, double s);
vec3 v3_add(vec3 a, vec3 b);
vec3 v3_unit(vec3 a);
vec3 m3_mul(const mat3 *A, vec3 x);
bool m3_inverse(const mat3 *A, mat3 *out);
quat q_mul(quat a, quat b);
quat q_conj(quat a);
quat q_axis_angle(vec3 axis, double angle);
vec3 q_rotate(quat q, vec3 v);   /* body -> inertial */
double q_error_angle(quat q, quat q_cmd);

/* Pick w_max from wheel momentum and alpha from wheel torque for this axis.
 * margin < 1 reserves capacity for gyroscopic compensation and the bias
 * momentum h0 stored in the wheels; brake_frac sets a = brake_frac * alpha so
 * the rate loop has acceleration headroom to track the braking curve.
 * rate_cap <= 0 means no operational cap. */
slew_gains slew_design(const spacecraft *sc, vec3 axis_body, double h0_norm,
                       double zeta, double wn, double margin, double brake_frac,
                       double rate_cap);

/* Rate command magnitude |w_cmd| = 2k L(e)/c for an error vector of norm en. */
double slew_rate_command(const slew_gains *g, double en);

/* Wheel torque on the body commanded by the given control law. */
vec3 slew_torque(ctrl_mode mode, const spacecraft *sc, const slew_gains *g,
                 const sc_state *s, quat q_cmd);

sim_result slew_simulate(ctrl_mode mode, const spacecraft *sc, const slew_gains *g,
                         sc_state s0, quat q_cmd, const sim_opts *o,
                         sim_logger log, void *ctx);

/* Eigenaxis angle history from the scalar equation
 *   th'' = sat_alpha( -2k s min(1, L(s)/s) - c th' ),  s = sin(th/2)
 * integrated with RK4 at step dt. Writes theta at multiples of out_dt. */
void slew_scalar_reference(const slew_gains *g, double th0, double t_end,
                           double dt, double out_dt, double *theta, int n_out);

/* Rest-to-rest bang-coast-bang time with the same alpha and w_max. */
double slew_bcb_time(double angle, double alpha, double w_max);

#endif
