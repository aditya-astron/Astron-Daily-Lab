/*
 * Astron Built this AP
 * Project 039 | Rocket propulsion | BATES solid motor internal ballistics
 * Author: Aditya (aditya-astron)
 * Date: 2026-11-15 IST
 *
 * Verification: isentropic-table values, hand-calculated c*, exact grain
 * geometry identities, closed-form equilibrium and blowdown solutions,
 * the linearised ignition time constant, mass conservation and the
 * vacuum-impulse identity, plus time-step convergence.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ballistics.h"
#include "motor_design.h"

static int fails = 0, checks = 0;

static void check(const char *name, double got, double want, double rtol)
{
    double err = fabs(got - want) / (fabs(want) > 0 ? fabs(want) : 1.0);
    int ok = err <= rtol;
    ++checks;
    if (!ok) ++fails;
    printf("[%s] %-46s got %-14.8g want %-14.8g rel %.2e\n",
           ok ? "PASS" : "FAIL", name, got, want, err);
}

static void check_true(const char *name, int cond)
{
    ++checks;
    if (!cond) ++fails;
    printf("[%s] %s\n", cond ? "PASS" : "FAIL", name);
}

static void test_isentropic(void)
{
    /* gamma 1.4, M 2: A/A* = (1.8/1.2)^3 / 2 = 1.6875 exactly; p/p0 = 1.8^-3.5
     * = 0.12780 (NACA 1135 tabulates 0.1278). */
    check("A/A* at M=2, gamma 1.4", area_ratio(2.0, 1.4), 1.6875, 1e-14);
    check("supersonic M from A/A* = 1.6875", mach_supersonic(1.6875, 1.4),
          2.0, 1e-12);
    check("pe/pc at M=2, gamma 1.4", pe_over_pc(1.6875, 1.4),
          pow(1.8, -3.5), 1e-12);
    check("pe/pc table value 0.1278", pe_over_pc(1.6875, 1.4), 0.1278, 5e-4);

    /* CF is maximised by the expansion ratio that matches pe to pa. */
    double g = 1.21, pc = 5e6, pa = PA_SEA_LEVEL;
    double eo = eps_for_optimum(pc, pa, g);
    check("eps_for_optimum gives pe == pa", pe_over_pc(eo, g) * pc, pa, 1e-10);
    double cf0 = thrust_coeff(eo, g, pc, pa);
    check_true("CF(eps_opt) > CF(0.9 eps_opt) and CF(1.1 eps_opt)",
               cf0 > thrust_coeff(0.9 * eo, g, pc, pa) &&
               cf0 > thrust_coeff(1.1 * eo, g, pc, pa));
    /* Vacuum CF rises monotonically toward its infinite-expansion limit;
     * for gamma 1.21 the approach is slow (pe/pc ~ eps^-gamma). */
    double cf_inf = sqrt(2 * g * g / (g - 1) * pow(2 / (g + 1), (g + 1) / (g - 1)));
    int mono = 1;
    double prev = 0;
    for (double e = 2; e < 1e13; e *= 10) {
        double c = thrust_coeff(e, g, 1.0, 0.0);
        if (c <= prev || c >= cf_inf) mono = 0;
        prev = c;
    }
    check_true("CF vacuum increases with eps and stays below its limit", mono);
    check("CF vacuum at eps 1e12 -> infinite-expansion limit",
          thrust_coeff(1e12, g, 1.0, 0.0), cf_inf, 1e-3);
}

static void test_cstar(void)
{
    Propellant p = design_propellant();
    /* Hand calculation: Gamma(1.21) = sqrt(1.21) (2/2.21)^(2.21/0.42)
     * = 1.1 * 0.59131 = 0.65044; sqrt(R Tc) = sqrt(346.436 * 3000)
     * = 1019.47 m/s; c* = 1567.4 m/s. */
    check("ideal c* vs hand calculation", cstar_ideal(&p), 1567.4, 1e-4);
    check("burn rate 7 mm/s at 6.9 MPa", p.a * pow(6.9e6, p.n), 0.007, 1e-12);
}

static double simpson(const Grain *g, double (*f)(const Grain *, double),
                      double a, double b, int n)
{
    double h = (b - a) / n, s = f(g, a) + f(g, b);
    for (int i = 1; i < n; ++i) s += (i % 2 ? 4.0 : 2.0) * f(g, a + i * h);
    return s * h / 3.0;
}

static void test_geometry(void)
{
    Grain g = design_grain();
    double w = grain_web(&g);
    /* Every bit of propellant is consumed by sweeping the burning surface
     * through the web: integral of Ab dx = V_prop(0). Ab is a quadratic in
     * x, so Simpson is exact up to round-off. */
    double integral = simpson(&g, grain_burn_area, 0.0, w * (1 - 1e-15), 2000);
    check("integral Ab dx == initial propellant volume", integral,
          grain_prop_volume(&g, 0.0), 1e-10);
    /* Ab = -dV/dx */
    double x = 0.37 * w, h = 1e-6;
    double dv = -(grain_prop_volume(&g, x + h) - grain_prop_volume(&g, x - h))
                / (2 * h);
    check("Ab == -dVprop/dx (central difference)", grain_burn_area(&g, x),
          dv, 1e-8);
    /* Neutral length equalises start and end area. */
    check("Ab(0) == Ab(web) at L = (3D+d)/2", grain_burn_area(&g, 0.0),
          grain_burn_area(&g, w * (1 - 1e-13)), 1e-9);
    /* Peak area sits at x* = (L - 2d)/6: dense scan. */
    double xbest = 0, Abest = 0;
    for (int i = 0; i <= 200000; ++i) {
        double xi = w * i / 200000.0 * (1 - 1e-12);
        double A = grain_burn_area(&g, xi);
        if (A > Abest) { Abest = A; xbest = xi; }
    }
    check("peak-area web location (L-2d)/6", xbest, (g.L - 2 * g.d) / 6.0,
          1e-4);
    /* Brute-force: among segment lengths, the neutral rule minimises the
     * max/min area ratio (for a concave quadratic the minimum is at an end,
     * and equal ends give the lowest range). */
    double Lbest = 0, rbest = 1e9;
    for (int i = 0; i <= 4000; ++i) {
        Grain t = g;
        t.L = 0.10 + 0.15 * i / 4000.0;
        double wt = grain_web(&t), Amin = 1e9, Amax = 0;
        for (int k = 0; k <= 400; ++k) {
            double A = grain_burn_area(&t, wt * k / 400.0 * (1 - 1e-12));
            if (A < Amin) Amin = A;
            if (A > Amax) Amax = A;
        }
        if (Amax / Amin < rbest) { rbest = Amax / Amin; Lbest = t.L; }
    }
    check("brute-force best L == (3D+d)/2", Lbest, g.L, 5e-4);
}

static void test_equilibrium_and_sizing(void)
{
    Propellant p = design_propellant();
    Grain g = design_grain();
    double At = size_throat(&p, &g, P_DESIGN_MAX);
    double xs = (g.L - 2 * g.d) / 6.0;
    double Kn = grain_burn_area(&g, xs) / At;
    check("equilibrium Pc at peak Kn == design MEOP",
          pc_equilibrium_exact(&p, Kn), P_DESIGN_MAX, 1e-10);
    /* The rho_g correction is small but real: factor (1 - rho_g/rho_p)^(1/(1-n)). */
    double rho_g = P_DESIGN_MAX / (p.R * p.Tc);
    check("rho_g correction factor", pc_equilibrium_exact(&p, Kn) /
          pc_equilibrium(&p, Kn), pow(1 - rho_g / p.rho, 1 / (1 - p.n)), 1e-9);

    Nozzle nz = { At, eps_for_optimum(0.95 * P_DESIGN_MAX, PA_SEA_LEVEL, p.gamma) };
    static TraceSample tr[8192];
    MotorResult r = simulate_motor(&p, &g, &nz, PA_SEA_LEVEL, 1e-5, P_IGNITER,
                                   tr, 8192, 100);
    check("simulated peak Pc == MEOP target", r.P_max, P_DESIGN_MAX, 2e-4);

    /* Mid-burn the chamber trails the moving quasi-steady equilibrium by
     * the fill time constant tau(t) = V / (R Tc ((1-n) At/c* + Ab r/(R Tc))),
     * i.e. Pc ~= Pc_eq - tau dPc_eq/dt. Late in the burn the area falls
     * fast and the chamber is ~5 L, so the raw lag reaches a few tenths of
     * a percent; the first-order lag model should remove nearly all of it. */
    double worst_raw = 0, worst_lag = 0;
    for (size_t i = 1; i + 1 < r.n_samples; ++i) {
        if (tr[i].t < 0.5 || tr[i + 1].t > r.t_burnout - 0.01) continue;
        double pe  = pc_equilibrium_exact(&p, tr[i].Ab / At);
        double pem = pc_equilibrium_exact(&p, tr[i - 1].Ab / At);
        double pep = pc_equilibrium_exact(&p, tr[i + 1].Ab / At);
        double dpe = (pep - pem) / (tr[i + 1].t - tr[i - 1].t);
        double V = grain_free_volume(&g, tr[i].x), RT = p.R * p.Tc;
        double rr = p.a * pow(tr[i].Pc, p.n);
        double tau = V / (RT * ((1 - p.n) * At / cstar(&p) + tr[i].Ab * rr / RT));
        double e1 = fabs(tr[i].Pc / pe - 1);
        double e2 = fabs(tr[i].Pc / (pe - tau * dpe) - 1);
        if (e1 > worst_raw) worst_raw = e1;
        if (e2 > worst_lag) worst_lag = e2;
    }
    printf("       mid-burn |Pc/Pc_eq - 1|: raw %.2e, lag-corrected %.2e\n",
           worst_raw, worst_lag);
    check_true("raw quasi-steady error below 1 %", worst_raw < 1e-2);
    check_true("lag-corrected equilibrium within 1e-4 (second-order residue)",
               worst_lag < 1e-4);

    /* Conservation: everything loaded leaves through the throat or is still
     * in the chamber. */
    check("mass: nozzle + chamber == loaded", r.m_nozzle + r.m_gas_end,
          r.m_prop, 1e-9);
    /* Vacuum impulse identity: integral Pc At dt = c* m_nozzle, CF constant. */
    check("I_vac == CF_vac c* m_nozzle", r.I_vac,
          thrust_coeff(nz.eps, p.gamma, 1, 0) * cstar(&p) * r.m_nozzle, 1e-9);
    /* Burnout time against web / mean burn rate from the trace. */
    check("free volume at burnout == case volume", r.V_free_burnout,
          0.25 * M_PI * g.D * g.D * g.Lcase, 1e-12);

    /* Step convergence: RK4 at 10 us vs 2.5 us. */
    MotorResult r2 = simulate_motor(&p, &g, &nz, PA_SEA_LEVEL, 2.5e-6,
                                    P_IGNITER, NULL, 0, 1);
    check("total impulse converged in dt", r.I_total, r2.I_total, 1e-6);
    check("burnout time converged in dt", r.t_burnout, r2.t_burnout, 1e-6);
}

static void test_blowdown(void)
{
    /* After web burnout the chamber is a fixed volume venting through a
     * choked throat: dPc/dt = -(R Tc At / (V c*)) Pc, so pure exponential. */
    Propellant p = design_propellant();
    Grain g = design_grain();
    Nozzle nz = { size_throat(&p, &g, P_DESIGN_MAX), 6.0 };
    static TraceSample tr[600000];
    MotorResult r = simulate_motor(&p, &g, &nz, PA_SEA_LEVEL, 1e-5, P_IGNITER,
                                   tr, 600000, 1);
    double k = p.R * p.Tc * nz.At / (r.V_free_burnout * cstar(&p));
    double worst = 0;
    int n = 0;
    for (size_t i = 0; i < r.n_samples; ++i) {
        if (tr[i].t <= r.t_burnout) continue;
        double want = r.P_burnout * exp(-k * (tr[i].t - r.t_burnout));
        double e = fabs(tr[i].Pc / want - 1);
        if (e > worst) worst = e;
        ++n;
    }
    printf("       tail-off: %d samples, tau = %.3f ms, worst rel err %.2e\n",
           n, 1e3 / k, worst);
    check_true("tail-off matches exp(-t/tau) to 1e-9", n > 100 && worst < 1e-9);
}

static void test_ignition_time_constant(void)
{
    /* Small-signal decay toward equilibrium. Linearising the chamber balance
     * at Pc_eq gives  d(dP)/dt = lambda dP  with
     *   lambda = (R Tc / V) [ -(1 - n) At / c* - Ab r / (R Tc) ].
     * Use L = 2d so dAb/dx = 0 at x = 0 and the equilibrium does not drift. */
    Propellant p = design_propellant();
    Grain g = { 4, 0.098, 0.035, 0.070, 4 * 0.070 + 0.01 };
    Nozzle nz = { size_throat(&p, &g, 4e6), 5.0 };
    double Ab = grain_burn_area(&g, 0), Peq = pc_equilibrium_exact(&p, Ab / nz.At);
    double V = grain_free_volume(&g, 0), RT = p.R * p.Tc;
    double r0 = p.a * pow(Peq, p.n);
    double lam = RT / V * (-(1 - p.n) * nz.At / cstar(&p) - Ab * r0 / RT);
    double tau = -1 / lam;
    double dt = 1e-6;
    static TraceSample tr[20000];
    double dP0 = 0.002 * Peq;
    simulate_motor(&p, &g, &nz, PA_SEA_LEVEL, dt, Peq + dP0, tr, 20000, 1);
    size_t i = (size_t)llround(tau / dt);
    double dP = tr[i].Pc - Peq;
    printf("       tau = %.4f ms, Pc_eq %.4f MPa\n", tau * 1e3, Peq / 1e6);
    check("deviation after one tau == e^-1", dP / dP0, exp(-1.0), 3e-3);
}

int main(void)
{
    test_isentropic();
    test_cstar();
    test_geometry();
    test_equilibrium_and_sizing();
    test_blowdown();
    test_ignition_time_constant();
    printf("\n%d/%d checks passed\n", checks - fails, checks);
    return fails ? EXIT_FAILURE : EXIT_SUCCESS;
}
