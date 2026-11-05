// Astron Built this AP
// Project 015 | Electronics and circuit analysis | Damped LC input filter vs. constant-power-load instability
// Author: Aditya (aditya-astron)
// Date: 2026-11-05 IST
#include <cmath>
#include <cstdio>

#include "filter.hpp"

static int failures = 0;

static void check(bool ok, const char* what, double got, double want) {
    std::printf("[%s] %-58s got %-12.6g want %.6g\n", ok ? "PASS" : "FAIL", what, got, want);
    if (!ok) ++failures;
}

static bool rel(double a, double b, double tol) { return std::abs(a - b) <= tol * std::abs(b); }

int main() {
    const double vin = 28.0;
    const double L = 10e-6, C = 22e-6, R0 = std::sqrt(L / C);
    const double f0 = 1.0 / (2.0 * M_PI * std::sqrt(L * C));

    // 1. Nodal solver reproduces the hand-derived Zout of the same network.
    {
        const FilterDesign d{L, C, 0.02, 0.4, 88e-6};
        const Netlist net = build_netlist(d);
        double worst = 0.0;
        for (double lf = 2.0; lf <= 6.0; lf += 0.01) {
            const double f = std::pow(10.0, lf);
            worst = std::max(worst, std::abs(net.impedance(1, f) - zout_closed_form(d, f)) /
                                        std::abs(zout_closed_form(d, f)));
        }
        check(worst < 1e-12, "nodal Zout vs closed form, max rel error", worst, 0.0);
    }

    // 2. Lightly damped LC: peak at f0 with |Zo| ~ Q*R0, Q = R0/RL.
    {
        const double RL = 0.005;
        const FilterDesign d{L, C, RL, 0.0, 0.0};
        double fpk;
        const double zpk = peak_zout(d, 100.0, 1e6, &fpk);
        check(rel(fpk, f0, 1e-3), "undamped resonance frequency f0 (Hz)", fpk, f0);
        check(rel(zpk, R0 * R0 / RL, 1e-3), "undamped peak |Zo| = R0^2/RL (Ohm)", zpk, R0 * R0 / RL);
    }

    // 3. Numerical minimax over Rd reproduces Erickson's closed-form optimum
    //    (lossless inductor) for several Cd/C ratios.
    for (double n : {1.0, 4.0, 10.0}) {
        FilterDesign d{L, C, 0.0, 0.0, n * C};
        const double rd = rd_minimax(d, 100.0, 1e6);
        d.Rd = rd;
        double fpk;
        const double zpk = peak_zout(d, 100.0, 1e6, &fpk);
        char buf[96];
        std::snprintf(buf, sizeof buf, "n=%-4g minimax Rd vs Erickson Rd_opt (Ohm)", n);
        check(rel(rd, rd_optimal(R0, n), 5e-3), buf, rd, rd_optimal(R0, n));
        std::snprintf(buf, sizeof buf, "n=%-4g optimal peak |Zo| (Ohm)", n);
        check(rel(zpk, zpeak_optimal(R0, n), 1e-4), buf, zpk, zpeak_optimal(R0, n));
        std::snprintf(buf, sizeof buf, "n=%-4g optimal peak frequency (Hz)", n);
        check(rel(fpk, fpeak_optimal(f0, n), 1e-2), buf, fpk, fpeak_optimal(f0, n));
    }

    // 4. Undamped filter + CPL: closed-form stability boundary.
    //    a1 = 0 gives Rn = L/(RL C); with v^2 - Vin v + RL P = 0 that yields
    //    v = Vin / (1 + RL^2 C / L) and P = v^2 RL C / L.
    {
        const double RL = 0.02;
        const FilterDesign d{L, C, RL, 0.0, 0.0};
        const double v = vin / (1.0 + RL * RL * C / L);
        const double p_crit = v * v * RL * C / L;
        check(rel(max_stable_power(d, vin, 5000.0), p_crit, 1e-6),
              "undamped Routh boundary vs closed form (W)", max_stable_power(d, vin, 5000.0), p_crit);
        const double below = simulate_growth(d, vin, 0.8 * p_crit, 0.05, 5e-3, 5e-9);
        const double above = simulate_growth(d, vin, 1.2 * p_crit, 0.05, 5e-3, 5e-9);
        check(below < 0.5, "sim at 0.8 P_crit decays (growth ratio)", below, 0.0);
        check(above > 2.0, "sim at 1.2 P_crit grows (growth ratio)", above, 0.0);
    }

    // 5. Damped design at 100 W: Middlebrook satisfied, Routh and simulation agree,
    //    and simulation brackets the Routh power limit.
    {
        const FilterDesign bare{L, C, 0.02, 0.0, 0.0};
        const FilterDesign d{L, C, 0.02, rd_optimal(R0, 4.0), 4.0 * C};
        const auto op = operating_point(d, vin, 100.0);
        const double zpk = peak_zout(d, 100.0, 1e6);
        check(op.rn / zpk > 2.0, "damped: Zin/peak Zout > 6 dB Middlebrook margin", op.rn / zpk, 2.0);
        check(!is_stable(bare, vin, 100.0), "undamped unstable at 100 W (Routh)", 0, 0);
        check(is_stable(d, vin, 100.0), "damped stable at 100 W (Routh)", 1, 1);
        check(std::isinf(simulate_growth(bare, vin, 100.0, 0.1, 5e-3, 5e-9)),
              "undamped at 100 W collapses in simulation", 1, 1);
        const double pmax = max_stable_power(d, vin, 5000.0);
        // Middlebrook is sufficient, not necessary: the true limit is above the
        // power where V^2/P equals the peak |Zo|. With the RL drop that power is
        // v = Vin / (1 + RL/Zpk), P = v^2 / Zpk.
        const double v_mb = vin / (1.0 + d.RL / zpk), p_mb = v_mb * v_mb / zpk;
        check(pmax > p_mb, "Routh P_max above Middlebrook-equal power (W)", pmax, p_mb);
        const double below = simulate_growth(d, vin, 0.9 * pmax, 0.05, 1e-2, 5e-9);
        const double above = simulate_growth(d, vin, 1.1 * pmax, 0.05, 1e-2, 5e-9);
        check(below < 1.0, "damped sim at 0.9 P_max decays", below, 0.0);
        check(above > 1.0, "damped sim at 1.1 P_max grows", above, 0.0);
    }

    std::printf("\n%s (%d failure%s)\n", failures ? "FAILED" : "ALL PASSED", failures,
                failures == 1 ? "" : "s");
    return failures ? 1 : 0;
}
