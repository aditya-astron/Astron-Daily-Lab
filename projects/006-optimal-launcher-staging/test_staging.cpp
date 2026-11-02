// Astron Built this AP
// Project 006 | Rocket propulsion and staging | Optimal mass staging of a multistage launcher
// Author: Aditya (aditya-astron)
// Date: 2026-11-02 IST
//
// Checks: hand-calculated single stage, closed form for identical stages,
// dv constraint and mass bookkeeping, a brute-force search that must not beat
// the Lagrange optimum, and infeasible-input handling.

#include <cmath>
#include <cstdio>
#include <vector>

#include "staging.hpp"

using namespace staging;

static int failures = 0;

static void check(bool ok, const char* what) {
    std::printf("[%s] %s\n", ok ? "PASS" : "FAIL", what);
    if (!ok) ++failures;
}

static bool close(double a, double b, double rel) {
    return std::fabs(a - b) <= rel * std::fabs(b);
}

int main() {
    // 1. Single stage by hand: c = 3000 m/s, eps = 0.1, dv = c ln 4 -> n = 4,
    //    m0/mPL = n(1-eps)/(1-eps n) = 4*0.9/0.6 = 6.
    {
        const std::vector<Stage> s = {{"A", 3000.0 / kG0, 0.1}};
        const Solution r = optimize(s, 3000.0 * std::log(4.0), 1000.0);
        check(r.ok && close(r.stages[0].n, 4.0, 1e-10), "single stage: n = 4");
        check(r.ok && close(r.gross, 6000.0, 1e-10), "single stage: m0 = 6000 kg");
        check(r.ok && close(r.stages[0].propellant, 4500.0, 1e-9) &&
                  close(r.stages[0].inert, 500.0, 1e-9),
              "single stage: 4500 kg propellant, 500 kg inert");
    }

    // 2. N identical stages: symmetry gives n = exp(dv / (N c)) and
    //    m0/mPL = [n(1-eps)/(1-eps n)]^N.
    {
        const double isp = 320.0, eps = 0.09, dv = 8500.0, mpl = 2500.0;
        for (int N = 2; N <= 4; ++N) {
            const std::vector<Stage> s(N, Stage{"S", isp, eps});
            const Solution r = optimize(s, dv, mpl);
            const double n = std::exp(dv / (N * isp * kG0));
            const double ratio = std::pow(n * (1 - eps) / (1 - eps * n), N);
            bool ok = r.ok && close(r.gross, mpl * ratio, 1e-9);
            for (int i = 0; ok && i < N; ++i) ok = close(r.stages[i].n, n, 1e-9);
            char msg[64];
            std::snprintf(msg, sizeof msg, "%d identical stages match closed form", N);
            check(ok, msg);
        }
    }

    const std::vector<Stage> stack = {
        {"S1", 300.0, 0.06}, {"S2", 345.0, 0.07}, {"S3", 450.0, 0.12}};
    const double dv = 9000.0, mpl = 10000.0;
    const Solution opt = optimize(stack, dv, mpl);
    check(opt.ok, "three-stage case solves");

    // 3. Constraint and bookkeeping: stage dv sums to target; masses rebuilt
    //    from inert + propellant reproduce ignition masses; rocket equation
    //    recomputed from masses gives the same stage dv.
    {
        double sum = 0.0, above = mpl, worst = 0.0;
        for (std::size_t k = stack.size(); k-- > 0;) {
            const auto& m = opt.stages[k];
            const double ign = above + m.inert + m.propellant;
            const double dvk = stack[k].c() * std::log(ign / (ign - m.propellant));
            worst = std::fmax(worst, std::fabs(ign - m.ignition) / ign);
            worst = std::fmax(worst, std::fabs(dvk - m.dv) / m.dv);
            sum += dvk;
            above = ign;
        }
        check(std::fabs(sum - dv) < 1e-6, "stage dv sums to 9000 m/s");
        check(worst < 1e-12, "mass bookkeeping and rocket equation consistent");
        check(close(above, opt.gross, 1e-12), "stacked masses equal lift-off mass");
    }

    // 4. Optimality condition 1/(1 - eps_i n_i) = eta c_i holds for every stage.
    {
        bool ok = true;
        for (std::size_t i = 0; i < stack.size(); ++i)
            ok = ok && close(1.0 / (1.0 - stack[i].eps * opt.stages[i].n),
                             opt.eta * stack[i].c(), 1e-10);
        check(ok, "first-order condition 1/(1-eps n) = eta c");
    }

    // 5. Brute force: sweep dv1, dv2 on a grid (dv3 from the constraint).
    //    No feasible point may be lighter than the optimum, and the best grid
    //    point should land within grid resolution of it.
    {
        double best = 1e300, best_dv1 = 0, best_dv2 = 0;
        const double step = 5.0;
        for (double d1 = step; d1 < dv; d1 += step)
            for (double d2 = step; d1 + d2 < dv; d2 += step) {
                const double d3 = dv - d1 - d2;
                const std::vector<double> n = {std::exp(d1 / stack[0].c()),
                                               std::exp(d2 / stack[1].c()),
                                               std::exp(d3 / stack[2].c())};
                const Solution r = evaluate(stack, n, mpl);
                if (r.ok && r.gross < best) {
                    best = r.gross;
                    best_dv1 = d1;
                    best_dv2 = d2;
                }
            }
        check(best >= opt.gross * (1.0 - 1e-12), "grid search never beats optimum");
        check(close(best, opt.gross, 1e-4), "grid minimum within 0.01% of optimum");
        check(std::fabs(best_dv1 - opt.stages[0].dv) <= 2 * step &&
                  std::fabs(best_dv2 - opt.stages[1].dv) <= 2 * step,
              "grid minimum dv split agrees with optimum");
        std::printf("       optimum %.1f kg, grid best %.1f kg\n", opt.gross, best);
    }

    // 6. Equal split is feasible here but heavier.
    {
        const Solution eq = equal_split(stack, dv, mpl);
        check(eq.ok && eq.gross > opt.gross, "equal dv split is heavier than optimum");
    }

    // 7. Infeasible inputs are reported, not silently solved.
    {
        check(!optimize(stack, max_delta_v(stack) + 1.0, mpl).ok,
              "dv above zero-payload limit rejected");
        const std::vector<Stage> one = {{"A", 300.0, 0.1}};
        check(!optimize(one, 7000.0, mpl).ok, "single stage to 7 km/s at eps 0.1 rejected");
        check(!evaluate(stack, {1.0, 2.0, 2.0}, mpl).ok, "n = 1 rejected by evaluate");
    }

    std::printf("\n%s (%d failure%s)\n", failures ? "FAILED" : "ALL PASSED", failures,
                failures == 1 ? "" : "s");
    return failures ? 1 : 0;
}
