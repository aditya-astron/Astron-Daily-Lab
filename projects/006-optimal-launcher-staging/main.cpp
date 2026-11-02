// Astron Built this AP
// Project 006 | Rocket propulsion and staging | Optimal mass staging of a multistage launcher
// Author: Aditya (aditya-astron)
// Date: 2026-11-02 IST
//
// Case study: size a three-stage launcher that puts 10 t into a 200 km LEO.
// Ideal budget: ~7.79 km/s circular speed plus ~1.6 km/s gravity, drag and
// steering losses, less ~0.4 km/s from an easterly launch at 28.5 deg
// latitude, rounded to a 9.0 km/s design dv. Isp values are burn-averaged
// (stage 1 sees part sea-level, part vacuum expansion).

#include <cstdio>
#include <vector>

#include "staging.hpp"

using staging::Solution;
using staging::Stage;

static void print(const char* title, const std::vector<Stage>& s,
                  const Solution& sol) {
    std::printf("\n%s\n", title);
    if (!sol.ok) {
        std::printf("  failed: %s\n", sol.error.c_str());
        return;
    }
    std::printf("  %-10s %6s %5s %7s %9s %11s %11s %11s\n", "stage", "Isp", "eps",
                "n", "dv [m/s]", "ign [kg]", "prop [kg]", "inert [kg]");
    double dv = 0.0;
    for (std::size_t i = 0; i < s.size(); ++i) {
        const auto& m = sol.stages[i];
        std::printf("  %-10s %6.0f %5.2f %7.4f %9.1f %11.0f %11.0f %11.0f\n",
                    s[i].name.c_str(), s[i].isp_s, s[i].eps, m.n, m.dv, m.ignition,
                    m.propellant, m.inert);
        dv += m.dv;
    }
    std::printf("  total dv %.1f m/s | lift-off %.0f kg | payload fraction %.4f\n",
                dv, sol.gross, sol.payload / sol.gross);
    if (sol.eta > 0.0) std::printf("  eta = %.6e s/m\n", sol.eta);
}

int main() {
    const double dv = 9000.0;
    const double payload = 10000.0;
    const std::vector<Stage> stack = {
        {"S1 RP-1", 300.0, 0.05},   // LOX/RP-1 booster, burn-averaged Isp
        {"S2 RP-1", 345.0, 0.07},   // vacuum LOX/RP-1
        {"S3 LH2", 450.0, 0.14},    // LOX/LH2, heavier tanks per kg propellant
    };

    std::printf("Optimal staging: payload %.0f kg, dv %.0f m/s (stack limit %.0f m/s)\n",
                payload, dv, staging::max_delta_v(stack));

    const Solution opt = staging::optimize(stack, dv, payload);
    print("Lagrange optimum", stack, opt);

    const Solution eq = staging::equal_split(stack, dv, payload);
    print("Equal dv split (baseline)", stack, eq);
    if (opt.ok && eq.ok)
        std::printf("\n  optimum saves %.0f kg of lift-off mass (%.1f%%)\n",
                    eq.gross - opt.gross, 100.0 * (eq.gross - opt.gross) / eq.gross);

    // Structural-ratio sensitivity: re-optimise with each eps raised by 0.01.
    std::printf("\nSensitivity, eps + 0.01 per stage (re-optimised):\n");
    for (std::size_t i = 0; i < stack.size(); ++i) {
        auto s = stack;
        s[i].eps += 0.01;
        const Solution r = staging::optimize(s, dv, payload);
        if (r.ok)
            std::printf("  %-10s lift-off %9.0f kg (%+.1f%%)\n", stack[i].name.c_str(),
                        r.gross, 100.0 * (r.gross - opt.gross) / opt.gross);
    }

    // How many stages? Identical stages, Isp 350 s, eps 0.08.
    std::printf("\nStage count, identical stages (Isp 350 s, eps 0.08):\n");
    for (int N = 1; N <= 6; ++N) {
        std::vector<Stage> s(N, Stage{"S", 350.0, 0.08});
        const Solution r = staging::optimize(s, dv, payload);
        if (r.ok)
            std::printf("  N=%d  lift-off %12.0f kg  payload fraction %.4f\n", N, r.gross,
                        payload / r.gross);
        else
            std::printf("  N=%d  infeasible (%s)\n", N, r.error.c_str());
    }
    return opt.ok ? 0 : 1;
}
