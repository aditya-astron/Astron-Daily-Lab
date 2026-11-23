// Astron Built this AP
// Project 058 | Rocket propulsion | Hydrazine blowdown system sizing with ullage-gas thermodynamics
// Author: Aditya (aditya-astron)
// Date: 2026-11-23 IST
//
// Study report: thruster operating line, blowdown-ratio trade, then the
// manoeuvre plan flown with the ullage energy balance for several wall heat
// transfer coefficients, and the loaded pressure that restores the floor.

#include <cstdio>

#include "design.hpp"

using namespace bd;

int main() {
    const Requirements req;
    const Thruster thr = baseline_thruster();
    const double rho_p = System{}.rho_p;
    const Budget bud = propellant_budget(req, thr, rho_p);

    std::printf("== Thruster (lumped feed + choked throat) ==\n");
    std::printf("throat d = %.3f mm, K = %.3e kg/s/sqrt(Pa), Isp = %.1f s\n",
                2e3 * std::sqrt(thr.at / kPi), thr.k, thr.isp());
    std::printf("%8s %8s %8s %10s\n", "Pf bar", "Pc bar", "F N", "mdot g/s");
    for (double pf : {5.5, 8.0, 11.0, 16.0, 22.0, 24.0})
        std::printf("%8.1f %8.3f %8.4f %10.4f\n", pf, thr.chamber_pressure(pf * kBar) / kBar,
                    thr.thrust(pf * kBar), 1e3 * thr.mdot(pf * kBar));

    std::printf("\n== Propellant budget ==\n");
    std::printf("m0 = %.1f kg, dV = %.0f m/s -> usable %.3f kg, loaded %.3f kg "
                "(%.0f %% residual), V_load = %.2f L\n",
                req.m_sc0, req.dv, bud.m_usable, bud.m_load, 100 * req.residual, 1e3 * bud.v_load);

    std::printf("\n== Blowdown-ratio optimum for pressure-sized tanks (min P0*Vt) ==\n");
    const double b1 = 1.0 + std::sqrt(bud.v_usable / bud.v_load);
    const double b1n = golden_min([&](double b) { return pv_metric(b, bud.v_load, bud.v_usable, 1.0); }, 1.05, 10.0);
    std::printf("isothermal, with residual: closed form B* = %.6f, golden search %.6f\n", b1, b1n);
    const double bg = golden_min([&](double b) { return pv_metric(b, bud.v_load, bud.v_usable, 1.4); }, 1.05, 10.0);
    std::printf("adiabatic (n = 1.4), with residual: golden search B* = %.4f "
                "(no-residual closed form %.4f)\n", bg, std::pow(2.4 / 1.4, 1.4));

    std::printf("\n== Isothermal sizing across blowdown ratio (Pend = %.1f bar) ==\n", req.p_min / kBar);
    std::printf("%6s %7s %7s %7s %9s %9s %9s %8s\n", "B", "P0 bar", "Vt L", "r mm", "t_p mm",
                "m_p kg", "m kg", "burn s");
    const TankMaterial mat;
    for (double b : {1.5, 1.985, 2.5, 3.0, 3.5, 4.0, 4.36}) {
        System s = isothermal_design(req, b, -1.0);
        ShellMass sm = shell_mass(s.vt, s.p0, mat);
        const double mp = mat.rho * 4 * kPi * sm.radius * sm.radius * sm.t_pressure;
        std::printf("%6.3f %7.2f %7.2f %7.1f %9.3f %9.3f %9.3f %8.0f%s\n", b, s.p0 / kBar,
                    1e3 * s.vt, 1e3 * sm.radius, 1e3 * sm.t_pressure, mp, sm.mass,
                    isothermal_burn_time(s), sm.gauge_limited ? "  gauge" : "");
    }
    std::printf("Every case is gauge-limited, so the P0*Vt optimum does not govern here:\n"
                "shell mass ~ Vt^(2/3) and the smallest tank (highest B the thrusters take) wins.\n");

    const double b_sel = 4.0;
    std::printf("\n== Manoeuvre plan with ullage thermodynamics, B = %.1f isothermal sizing ==\n", b_sel);
    std::printf("%.0f s burns + %.0f s coasts, 4 thrusters firing together\n",
                req.burn, req.coast);
    std::printf("%10s %6s %10s %10s %10s %10s\n", "h W/m2K", "burns", "Pmin bar", "Tmin K",
                "fire s", "dV m/s");
    for (double h : {-1.0, 50.0, 5.0, 1.0, 0.0}) {
        System s = isothermal_design(req, b_sel, h);
        MissionResult m = fly_mission(s, req.burn, req.coast);
        if (h < 0)
            std::printf("%10s", "isotherm");
        else
            std::printf("%10.1f", h);
        std::printf(" %6zu %10.3f %10.1f %10.0f %10.2f\n", m.burns.size(), m.p_min / kBar,
                    m.t_gas_min, m.fire_time, m.final_state.dv);
    }

    System s5 = isothermal_design(req, b_sel, 5.0);
    std::printf("\ntau = mg cv / (h A_g) at start (h = 5): %.0f s\n",
                s5.gas_mass() * s5.gas.cv() / (s5.h * s5.wetted_area(s5.vg0)));
    MissionResult m5 = fly_mission(s5, req.burn, req.coast);
    // Isothermal pressure at the same burn boundaries, for comparison.
    MissionResult mi = fly_mission(isothermal_design(req, b_sel, -1.0), req.burn, req.coast);
    std::printf("\nPer-burn history, h = 5 W/m2K\n");
    std::printf("%4s %8s %8s %9s %9s %8s %9s %8s\n", "#", "t0 h", "dur s", "Pst bar", "Pend bar",
                "Tend K", "Piso bar", "dV m/s");
    for (size_t i = 0; i < m5.burns.size(); ++i) {
        const BurnRecord& b = m5.burns[i];
        const double piso = i < mi.burns.size() ? mi.burns[i].p_end / kBar : 0.0;
        std::printf("%4zu %8.2f %8.1f %9.3f %9.3f %8.1f %9.3f %8.2f\n", i + 1, b.t_start / 3600,
                    b.duration, b.p_start / kBar, b.p_end / kBar, b.t_gas_end, piso, b.dv);
    }

    std::printf("\n== Loaded pressure that holds the %.1f bar floor (same tank) ==\n", req.p_min / kBar);
    std::printf("%10s %10s %10s %12s\n", "h W/m2K", "P0 bar", "dP0 bar", "eff. B");
    for (double h : {50.0, 5.0, 1.0, 0.0}) {
        System s = isothermal_design(req, b_sel, h);
        const double p0 = resize_p0(s, req);
        if (p0 < 0)
            std::printf("%10.1f %10s  exceeds %.0f bar thruster limit\n", h, "--", req.p_max / kBar);
        else
            std::printf("%10.1f %10.3f %10.3f %12.3f\n", h, p0 / kBar, (p0 - s.p0) / kBar,
                        p0 / req.p_min);
    }

    std::printf("\n== Flow-rate sweep: when the ullage stops being isothermal ==\n");
    std::printf("Same tank and load, h = 5 W/m2K, one continuous burn to depletion. The gas\n"
                "lags the wall by about tau * (adiabatic cooling rate), so what matters is the\n"
                "expulsion time against tau = %.0f s, not the burn plan.\n",
                s5.gas_mass() * s5.gas.cv() / (s5.h * s5.wetted_area(s5.vg0)));
    std::printf("%8s %10s %10s %10s %10s %10s %11s\n", "F/thr N", "t_exp s", "Pmin bar",
                "Tmin K", "P0 bar", "Vg0 L", "Vt L @22bar");
    for (double f : {1.0, 5.0, 10.0, 20.0, 50.0}) {
        Requirements r = req;
        r.thrust_nom = f;
        r.burn = 1e9;
        System s = isothermal_design(r, b_sel, 5.0);
        MissionResult m = fly_mission(s, r.burn, r.coast);
        const double p0 = resize_p0(s, r);
        const double vg0 = resize_vg0(s, r);
        std::printf("%8.0f %10.0f %10.3f %10.1f ", f, m.fire_time, m.p_min / kBar, m.t_gas_min);
        if (p0 < 0) std::printf("%10s ", "> 24");
        else std::printf("%10.3f ", p0 / kBar);
        std::printf("%10.3f %11.2f\n", 1e3 * vg0, 1e3 * (s.vt - s.vg0 + vg0));
    }
    return 0;
}
