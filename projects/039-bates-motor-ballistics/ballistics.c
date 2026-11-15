/*
 * Astron Built this AP
 * Project 039 | Rocket propulsion | BATES solid motor internal ballistics
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-15 IST
 */
#include "ballistics.h"

#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---------------------------------------------------------------- gas -- */

/* Vandenkerckhove function: mdot = Gamma * Pc * At / sqrt(R Tc). */
static double vdk(double g)
{
    return sqrt(g) * pow(2.0 / (g + 1.0), (g + 1.0) / (2.0 * (g - 1.0)));
}

double cstar_ideal(const Propellant *p)
{
    return sqrt(p->R * p->Tc) / vdk(p->gamma);
}

double cstar(const Propellant *p) { return p->eta_c * cstar_ideal(p); }

double area_ratio(double M, double g)
{
    double t = (2.0 / (g + 1.0)) * (1.0 + 0.5 * (g - 1.0) * M * M);
    return pow(t, (g + 1.0) / (2.0 * (g - 1.0))) / M;
}

/* A/A* is monotonic on the supersonic branch, so bisection in M is safe
 * and needs no derivative or starting guess. */
double mach_supersonic(double eps, double g)
{
    double lo = 1.0, hi = 2.0;
    while (area_ratio(hi, g) < eps) hi *= 2.0;
    for (int i = 0; i < 200; ++i) {
        double mid = 0.5 * (lo + hi);
        if (area_ratio(mid, g) < eps) lo = mid; else hi = mid;
        if (hi - lo < 1e-14 * hi) break;
    }
    return 0.5 * (lo + hi);
}

double pe_over_pc(double eps, double g)
{
    double M = mach_supersonic(eps, g);
    return pow(1.0 + 0.5 * (g - 1.0) * M * M, -g / (g - 1.0));
}

/* Ideal CF = momentum term + pressure term (Sutton eq. 3-30). */
/* pe/pc is fixed by geometry, so the hot loop passes it in precomputed. */
static double thrust_coeff_pr(double pr, double eps, double g, double pc,
                              double pa)
{
    double k  = 2.0 * g * g / (g - 1.0)
              * pow(2.0 / (g + 1.0), (g + 1.0) / (g - 1.0));
    double mom = sqrt(k * (1.0 - pow(pr, (g - 1.0) / g)));
    return mom + (pr - pa / pc) * eps;
}

double thrust_coeff(double eps, double g, double pc, double pa)
{
    return thrust_coeff_pr(pe_over_pc(eps, g), eps, g, pc, pa);
}

/* -------------------------------------------------------------- grain -- */

double grain_web(const Grain *g)
{
    /* Burnout is whichever comes first: core reaching the case wall or the
     * two end faces meeting. */
    double w_radial = 0.5 * (g->D - g->d);
    double w_axial  = 0.5 * g->L;
    return w_radial < w_axial ? w_radial : w_axial;
}

double grain_burn_area(const Grain *g, double x)
{
    if (x >= grain_web(g)) return 0.0;
    double dc = g->d + 2.0 * x;
    double Ls = g->L - 2.0 * x;
    double core = M_PI * dc * Ls;
    double ends = 0.5 * M_PI * (g->D * g->D - dc * dc);
    return g->N * (core + ends);
}

double grain_prop_volume(const Grain *g, double x)
{
    if (x >= grain_web(g)) return 0.0;
    double dc = g->d + 2.0 * x;
    double Ls = g->L - 2.0 * x;
    return g->N * 0.25 * M_PI * (g->D * g->D - dc * dc) * Ls;
}

double grain_free_volume(const Grain *g, double x)
{
    return 0.25 * M_PI * g->D * g->D * g->Lcase - grain_prop_volume(g, x);
}

double grain_neutral_length(double D, double d) { return 0.5 * (3.0 * D + d); }

/* ------------------------------------------------------- equilibrium -- */

double pc_equilibrium(const Propellant *p, double Kn)
{
    return pow(p->a * p->rho * cstar(p) * Kn, 1.0 / (1.0 - p->n));
}

/* Steady state of the chamber balance: gas generated = Ab r rho_p, but the
 * cavity the burning surface leaves behind fills with gas at rho_g, so only
 * Ab r (rho_p - rho_g) leaves through the throat. */
double pc_equilibrium_exact(const Propellant *p, double Kn)
{
    double P = pc_equilibrium(p, Kn);
    for (int i = 0; i < 100; ++i) {
        double rho_g = P / (p->R * p->Tc);
        double Pn = pow(p->a * cstar(p) * Kn * (p->rho - rho_g),
                        1.0 / (1.0 - p->n));
        if (fabs(Pn - P) < 1e-12 * P) { P = Pn; break; }
        P = Pn;
    }
    return P;
}

static double grain_max_area(const Grain *g)
{
    /* Ab(x) is a concave quadratic: dAb/dx = N pi (2L - 4d - 12x). */
    double w = grain_web(g);
    double best = grain_burn_area(g, 0.0);
    double xs = (g->L - 2.0 * g->d) / 6.0;
    double cand[2] = { xs, w * (1.0 - 1e-12) };
    for (int i = 0; i < 2; ++i)
        if (cand[i] > 0.0 && cand[i] < w) {
            double A = grain_burn_area(g, cand[i]);
            if (A > best) best = A;
        }
    return best;
}

double size_throat(const Propellant *p, const Grain *g, double Pmax)
{
    double rho_g = Pmax / (p->R * p->Tc);
    double Kn = pow(Pmax, 1.0 - p->n) / (p->a * cstar(p) * (p->rho - rho_g));
    return grain_max_area(g) / Kn;
}

/* --------------------------------------------------------- simulation -- */

enum { S_X, S_PC, S_MNOZ, S_I, S_IVAC, S_DIM };

typedef struct {
    const Propellant *p;
    const Grain *g;
    const Nozzle *nz;
    double cs, RT, cf_vac, pr, w, pa;
    int burnt;
} Ctx;

/* Chamber gas mass m = Pc V / (R Tc) with V the free volume. Differentiating,
 *   V/(R Tc) dPc/dt + rho_g dV/dt = rho_p Ab r - Pc At / c*,   dV/dt = Ab r
 * so dPc/dt = (R Tc / V) [ Ab r (rho_p - rho_g) - Pc At / c* ]. */
static void deriv(const Ctx *c, const double *s, double *ds)
{
    double Pc = s[S_PC] > 0.0 ? s[S_PC] : 0.0;
    double x  = c->burnt ? c->w : s[S_X];
    double Ab = c->burnt ? 0.0 : grain_burn_area(c->g, x);
    double r  = c->burnt ? 0.0 : c->p->a * pow(Pc, c->p->n);
    double V  = grain_free_volume(c->g, x);
    double rho_g = Pc / c->RT;
    double mdot  = Pc * c->nz->At / c->cs;

    ds[S_X]    = r;
    ds[S_PC]   = c->RT / V * (Ab * r * (c->p->rho - rho_g) - mdot);
    ds[S_MNOZ] = mdot;
    double F = Pc > 0.0
        ? thrust_coeff_pr(c->pr, c->nz->eps, c->p->gamma, Pc, c->pa)
              * Pc * c->nz->At
        : 0.0;
    ds[S_I]    = F > 0.0 ? F : 0.0;   /* no reverse thrust when over-expanded */
    ds[S_IVAC] = c->cf_vac * Pc * c->nz->At;
}

static void rk4(const Ctx *c, const double *s, double h, double *out)
{
    double k1[S_DIM], k2[S_DIM], k3[S_DIM], k4[S_DIM], tmp[S_DIM];
    deriv(c, s, k1);
    for (int i = 0; i < S_DIM; ++i) tmp[i] = s[i] + 0.5 * h * k1[i];
    deriv(c, tmp, k2);
    for (int i = 0; i < S_DIM; ++i) tmp[i] = s[i] + 0.5 * h * k2[i];
    deriv(c, tmp, k3);
    for (int i = 0; i < S_DIM; ++i) tmp[i] = s[i] + h * k3[i];
    deriv(c, tmp, k4);
    for (int i = 0; i < S_DIM; ++i)
        out[i] = s[i] + h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
}

static double thrust_now(const Ctx *c, double Pc)
{
    if (Pc <= 0.0) return 0.0;
    double F = thrust_coeff_pr(c->pr, c->nz->eps, c->p->gamma, Pc, c->pa)
             * Pc * c->nz->At;
    return F > 0.0 ? F : 0.0;
}

MotorResult simulate_motor(const Propellant *p, const Grain *g,
                           const Nozzle *nz, double pa, double dt,
                           double P_ign, TraceSample *trace,
                           size_t trace_cap, size_t trace_every)
{
    Ctx c = { p, g, nz, cstar(p), p->R * p->Tc,
              thrust_coeff(nz->eps, p->gamma, 1.0, 0.0),
              pe_over_pc(nz->eps, p->gamma), grain_web(g), pa, 0 };
    MotorResult res = {0};
    res.m_prop = p->rho * grain_prop_volume(g, 0.0);

    double s[S_DIM] = { 0.0, P_ign, 0.0, 0.0, 0.0 };
    double m_gas0 = P_ign * grain_free_volume(g, 0.0) / c.RT;
    double t = 0.0;
    size_t step = 0, n = 0;
    /* Below about twice ambient the nozzle unchokes and the c* mass-flow law
     * no longer holds, so the run stops there (or at 0.1 % of peak in vacuum). */
    double P_floor = pa > 0.0 ? 2.0 * pa : 0.0;

    for (;;) {
        if (trace && n < trace_cap && step % trace_every == 0) {
            double x = c.burnt ? c.w : s[S_X];
            trace[n++] = (TraceSample){ t, x, s[S_PC], thrust_now(&c, s[S_PC]),
                                        c.burnt ? 0.0 : grain_burn_area(g, x) };
        }
        double nxt[S_DIM];
        rk4(&c, s, dt, nxt);

        if (!c.burnt && nxt[S_X] >= c.w) {
            /* Locate web burnout inside the step by bisection on the RK4
             * step length, so the area discontinuity is never straddled. */
            double lo = 0.0, hi = dt;
            for (int i = 0; i < 60; ++i) {
                double mid = 0.5 * (lo + hi);
                rk4(&c, s, mid, nxt);
                if (nxt[S_X] < c.w) lo = mid; else hi = mid;
            }
            rk4(&c, s, lo, nxt);
            for (int i = 0; i < S_DIM; ++i) s[i] = nxt[i];
            t += lo;
            s[S_X] = c.w;
            c.burnt = 1;
            res.t_burnout = t;
            res.P_burnout = s[S_PC];
            res.V_free_burnout = grain_free_volume(g, c.w);
            ++step;
            continue;
        }

        for (int i = 0; i < S_DIM; ++i) s[i] = nxt[i];
        t += dt;
        ++step;
        if (s[S_PC] > res.P_max) res.P_max = s[S_PC];
        double F = thrust_now(&c, s[S_PC]);
        if (F > res.F_max) res.F_max = F;

        double floor_now = P_floor > 0.0 ? P_floor : 1e-3 * res.P_max;
        if (c.burnt && s[S_PC] < floor_now) break;
        if (t > 1e3) break;   /* guard against a non-igniting design */
    }

    res.t_end    = t;
    res.I_total  = s[S_I];
    res.I_vac    = s[S_IVAC];
    res.m_nozzle = s[S_MNOZ];
    /* Report gas left relative to the initial igniter gas so that
     * m_nozzle + m_gas_end == m_prop exactly for a conservative scheme. */
    res.m_gas_end = s[S_PC] * grain_free_volume(g, c.w) / c.RT - m_gas0;
    res.n_samples = n;
    return res;
}
