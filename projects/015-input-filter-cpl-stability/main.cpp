// Astron Built this AP
// Project 015 | Electronics and circuit analysis | Damped LC input filter vs. constant-power-load instability
// Author: Aditya (aditya-astron)
// Date: 2026-11-05 IST
//
// A 100 W point-of-load converter on a 28 V spacecraft bus needs an LC input
// filter for conducted-EMI limits. A regulated converter draws constant power,
// so its small-signal input impedance is a negative resistance -V^2/P, and an
// undamped LC in front of it becomes an oscillator. This sizes an Rd-Cd
// damping leg and checks the result three ways: Middlebrook impedance margin,
// Routh-Hurwitz on the linearised system, and a nonlinear time simulation.
#include <cmath>
#include <cstdio>

#include "filter.hpp"

int main() {
    const double vin = 28.0, p = 100.0;
    const double L = 10e-6, C = 22e-6, RL = 0.020, n = 4.0;
    const double f_lo = 100.0, f_hi = 1e6;

    const FilterDesign bare{L, C, RL, 0.0, 0.0};
    const double R0 = r0(bare);
    const double f0 = 1.0 / (2.0 * M_PI * std::sqrt(L * C));
    const FilterDesign damped{L, C, RL, rd_optimal(R0, n), n * C};

    const auto op = operating_point(bare, vin, p);
    std::printf("Bus %.0f V, converter %.0f W, L = %.0f uH, C = %.0f uF, RL = %.0f mOhm\n",
                vin, p, L * 1e6, C * 1e6, RL * 1e3);
    std::printf("R0 = sqrt(L/C) = %.3f Ohm, f0 = %.2f kHz\n", R0, f0 / 1e3);
    std::printf("Converter input: v = %.3f V, Zin = -V^2/P = -%.2f Ohm\n\n", op.v, op.rn);

    std::printf("Damping leg (Erickson optimum, n = Cd/C = %.0f): Cd = %.0f uF, Rd = %.4f Ohm\n",
                n, damped.Cd * 1e6, damped.Rd);
    std::printf("Numerical minimax Rd (with RL included):         Rd = %.4f Ohm\n\n",
                rd_minimax(damped, f_lo, f_hi));

    std::printf("%-10s %12s %12s %16s %8s %14s %16s\n", "filter", "peak|Zo| Ohm",
                "f_peak kHz", "Middlebrook dB", "Routh", "P_max W", "sim growth x");
    for (const auto* d : {&bare, &damped}) {
        double fpk;
        const double zpk = peak_zout(*d, f_lo, f_hi, &fpk);
        const double margin = 20.0 * std::log10(op.rn / zpk);
        const bool st = is_stable(*d, vin, p);
        const double pmax = max_stable_power(*d, vin, 5000.0);
        const double g = simulate_growth(*d, vin, p, 0.1, 5e-3, 5e-9);
        std::printf("%-10s %12.3f %12.2f %16.1f %8s %14.1f ", d == &bare ? "undamped" : "damped",
                    zpk, fpk / 1e3, margin, st ? "stable" : "UNSTABLE", pmax);
        if (std::isinf(g)) std::printf("%16s\n", "bus collapse");
        else std::printf("%16.3g\n", g);
    }

    std::printf("\n|Zout| sweep (Ohm) against |Zin| = %.2f Ohm\n", op.rn);
    std::printf("%10s %12s %12s\n", "f (kHz)", "undamped", "damped");
    const Netlist nb = build_netlist(bare), nd = build_netlist(damped);
    for (double f : {1e3, 3e3, 6e3, 8e3, 10e3, 10.73e3, 12e3, 15e3, 30e3, 100e3}) {
        std::printf("%10.2f %12.3f %12.3f\n", f / 1e3, std::abs(nb.impedance(1, f)),
                    std::abs(nd.impedance(1, f)));
    }
    return 0;
}
