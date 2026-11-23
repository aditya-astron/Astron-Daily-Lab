// Astron Built this AP
// Project 058 | Rocket propulsion | Hydrazine blowdown system sizing with ullage-gas thermodynamics
// Author: Aditya (aditya-astron)
// Date: 2026-11-23 IST
//
// Verification against closed forms and hand numbers. Every check prints its
// error; the binary exits non-zero if any fails.

#include <cstdio>

#include "design.hpp"

using namespace bd;

static int failures = 0;

static void check(const char* name, double got, double want, double rel_tol) {
    const double err = std::fabs(got - want) / std::max(std::fabs(want), 1e-300);
    const bool ok = err <= rel_tol;
    if (!ok) ++failures;
    std::printf("[%s] %-58s got %.10g want %.10g rel %.2e\n", ok ? "PASS" : "FAIL", name, got,
                want, err);
}

int main() {
    const Requirements req;
    const Thruster thr = baseline_thruster();

    // 1. Thruster operating line.
    check("nominal thrust at 22 bar", thr.thrust(22 * kBar), 1.0, 1e-12);
    check("nominal chamber pressure at 22 bar", thr.chamber_pressure(22 * kBar), 11 * kBar, 1e-12);
    check("Isp = c* Cf / g0 (hand: 223.44 s)", thr.isp(), 223.44, 1e-4);
    {
        // Independent root: bisection on the flow balance K sqrt(Pf-Pc) = Pc At/c*.
        const double pf = 7.3 * kBar;
        double lo = 0, hi = pf;
        for (int i = 0; i < 200; ++i) {
            const double pc = 0.5 * (lo + hi);
            if (thr.k * std::sqrt(pf - pc) > pc * thr.at / thr.cstar) lo = pc; else hi = pc;
        }
        check("Pc quadratic root vs bisection at 7.3 bar", thr.chamber_pressure(pf), 0.5 * (lo + hi), 1e-12);
        // Asymptotes: low feed pressure -> throat-limited Pc -> Pf; high -> Pc ~ sqrt(Pf).
        const double pf_lo = 1.0;  // Pa
        check("low-Pf asymptote Pc/Pf -> 1", thr.chamber_pressure(pf_lo) / pf_lo, 1.0, 1e-5);
        const double pf_hi = 1e12;
        check("high-Pf asymptote Pc -> K c* sqrt(Pf)/At", thr.chamber_pressure(pf_hi),
              thr.k * thr.cstar * std::sqrt(pf_hi) / thr.at, 1e-3);
    }

    // 2. Propellant budget by hand: 180 (1 - exp(-120/(223.44*9.80665))) = 9.5925 kg.
    const Budget bud = propellant_budget(req, thr, 1004.0);
    check("usable propellant (hand 9.5925 kg)", bud.m_usable, 9.5925, 2e-4);

    // 3. Blowdown-ratio optima.
    check("n=1 with residual: B* = 1 + sqrt(Vu/Vl)",
          golden_min([&](double b) { return pv_metric(b, bud.v_load, bud.v_usable, 1.0); }, 1.05, 10.0),
          1.0 + std::sqrt(bud.v_usable / bud.v_load), 1e-7);
    for (double n : {1.0, 1.2, 1.4, 1.67}) {
        char name[64];
        std::snprintf(name, sizeof name, "no residual, n=%.2f: B* = ((n+1)/n)^n", n);
        check(name, golden_min([&](double b) { return pv_metric(b, 1.0, 1.0, n); }, 1.05, 10.0),
              std::pow((n + 1) / n, n), 1e-7);
    }

    // 4. Shell sizing by hand: 1.5 * 22e5 * 0.1460 / (2 * 895e6) = 0.2692 mm.
    {
        System s = isothermal_design(req, 4.0, -1.0);
        ShellMass sm = shell_mass(s.vt, s.p0, TankMaterial{});
        check("tank radius from 13.03 L (hand 146.0 mm)", sm.radius, 0.14602, 2e-4);
        check("membrane thickness (hand 0.2692 mm)", sm.t_pressure, 1.5 * 22e5 * sm.radius / (2 * 895e6), 1e-12);
        check("gauge-limited mass 4 pi r^2 t rho", sm.mass,
              4430.0 * 4 * kPi * sm.radius * sm.radius * 0.6e-3, 1e-12);
    }

    // 5. Isothermal continuous burn: ODE depletion time vs adaptive quadrature,
    //    and the depletion event lands on Vg_end.
    {
        System s = isothermal_design(req, 4.0, -1.0);
        MissionResult m = fly_mission(s, 1e9, 0.0);
        check("isothermal burn time: ODE vs quadrature", m.fire_time, isothermal_burn_time(s), 1e-9);
        check("depletion event at Vg0 + Vu", m.final_state.vg, s.vg0 + s.v_usable, 1e-12);
        check("isothermal end pressure = P0/B", pressure(s, m.final_state), req.p_min, 1e-12);
        check("expelled mass = rho dVg", s.m_sc0 - m.final_state.m_sc,
              s.rho_p * (m.final_state.vg - s.vg0), 1e-12);
        check("integrated dV vs Isp g0 ln(m0/mf)", m.final_state.dv,
              thr.isp() * kG0 * std::log(s.m_sc0 / m.final_state.m_sc), 1e-10);
        check("dV at depletion = requirement", m.final_state.dv, req.dv, 1e-10);
    }

    // 6. Adiabatic limit (h = 0): Poisson invariants along the whole plan.
    {
        System s = isothermal_design(req, 4.0, 0.0);
        MissionResult m = fly_mission(s, req.burn, req.coast);
        const State& x = m.final_state;
        const double g = s.gas.gamma;
        check("adiabatic P Vg^gamma invariant", pressure(s, x) * std::pow(x.vg, g),
              s.p0 * std::pow(s.vg0, g), 1e-10);
        check("adiabatic T Vg^(gamma-1) invariant", x.temp * std::pow(x.vg, g - 1),
              s.tw * std::pow(s.vg0, g - 1), 1e-10);
        check("adiabatic floor = 22 bar / 4^1.4 (hand 3.1590 bar)", m.p_min, 22 * kBar / std::pow(4.0, 1.4), 1e-10);
    }

    // 7. Coast relaxation at fixed volume: T = Tw + (T0 - Tw) exp(-t/tau).
    {
        System s = isothermal_design(req, 4.0, 5.0);
        State x = initial_state(s);
        x.vg = s.vg0 + 0.5 * s.v_usable;
        x.temp = 240.0;
        const double tau = s.gas_mass() * s.gas.cv() / (s.h * s.wetted_area(x.vg));
        const double t_end = 150.0;
        for (int i = 0; i < 600; ++i) x = rk4_step(s, x, 0.25, false);
        check("coast recovery vs exponential, 150 s", x.temp,
              s.tw + (240.0 - s.tw) * std::exp(-t_end / tau), 1e-9);
    }

    // 8. Quasi-steady thermal lag: with t_exp >> tau the gas trails the wall by
    //    dT ~ P Qv / (h A_g) (expansion work balanced by wall heating). At
    //    depletion of a slow 1 N continuous burn the ODE should agree to ~1 %.
    {
        System s = isothermal_design(req, 4.0, 5.0);
        MissionResult m = fly_mission(s, 1e9, 0.0);
        const State& x = m.final_state;
        const double p = pressure(s, x);
        const double qv = s.n_thrusters * s.thr.mdot(p) / s.rho_p;
        const double dt_qs = p * qv / (s.h * s.wetted_area(x.vg));
        check("quasi-steady lag Tw - T at depletion", s.tw - x.temp, dt_qs, 1e-2);
    }

    // 9. Very large h reproduces the isothermal limit (tau ~ 0.1 s, dt below it).
    {
        Requirements r = req;
        System si = isothermal_design(r, 4.0, -1.0);
        System sh = isothermal_design(r, 4.0, 2000.0);
        MissionResult mi = fly_mission(si, r.burn, r.coast, 0.02);
        MissionResult mh = fly_mission(sh, r.burn, r.coast, 0.02);
        check("h -> inf converges to isothermal fire time", mh.fire_time, mi.fire_time, 1e-4);
    }

    // 10. Re-sized loaded pressure holds the floor; extra ullage does the same.
    {
        System s = isothermal_design(req, 4.0, 1.0);
        const double p0 = resize_p0(s, req);
        System s2 = s;
        s2.p0 = p0;
        check("resized P0 floor = 5.5 bar", fly_mission(s2, req.burn, req.coast).p_min, req.p_min, 1e-9);
        Requirements r = req;
        r.thrust_nom = 50.0;
        r.burn = 1e9;
        System s3 = isothermal_design(r, 4.0, 5.0);
        const double vl = s3.vt - s3.vg0;
        const double vg0 = resize_vg0(s3, r);
        s3.vg0 = vg0;
        s3.vt = vl + vg0;
        check("resized ullage floor = 5.5 bar (50 N case)", fly_mission(s3, r.burn, r.coast).p_min,
              req.p_min, 1e-9);
    }

    std::printf("\n%s: %d failure(s)\n", failures ? "FAILED" : "ALL PASSED", failures);
    return failures ? 1 : 0;
}
