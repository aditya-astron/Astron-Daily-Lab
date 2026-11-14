// Astron Built this AP
// Project 037 | ECLSS | Water recovery sorbent bed breakthrough and LUB scale-up
// Author: Aditya (aditya-astron)
// Date: 2026-11-14 IST
//
// Verification of the characteristic-grid solver against closed forms.

#include <cmath>
#include <cstdio>

#include "sorbent_bed.hpp"

static int failures = 0;

static void check(const char* what, double got, double want, double tol) {
    const double err = std::fabs(got - want);
    const bool ok = err <= tol;
    if (!ok) ++failures;
    std::printf("[%s] %-58s got %.6g want %.6g |err| %.2e tol %.1e\n", ok ? "PASS" : "FAIL", what, got,
                want, err, tol);
}

// Linear-isotherm bed with xi = 20 transfer units.
static bed::Bed linear_bed() {
    bed::Bed B{};
    B.length = 10.0;
    B.eps = 0.4;
    B.rho_b = 0.45;
    B.u = 20.0;
    B.k = 0.5;
    B.c0 = 0.02;
    // K chosen so xi = rho_b k K L / (eps u) = 20
    const double K = 20.0 * B.eps * B.u / (B.rho_b * B.k * B.length);
    B.iso = {K, 1.0, true};
    return B;
}

static double max_err_vs_thomas(const bed::Bed& B, int nz, int ntau, double tau_end) {
    const double xi = B.rho_b * B.k * B.iso.qm * B.iso.b * B.length / (B.eps * B.u);
    const auto bt = bed::simulate(B, nz, tau_end, ntau);
    double e = 0.0;
    for (size_t j = 0; j < bt.x.size(); j += 4)
        e = std::fmax(e, std::fabs(bt.x[j] - bed::thomas_J(xi, B.k * bt.tau[j])));
    return e;
}

int main() {
    std::printf("-- Anzelius/Thomas J function --\n");
    check("J(xi, 0) = exp(-xi), xi = 7", bed::thomas_J(7.0, 0.0), std::exp(-7.0), 1e-13);
    // J(x,y) + J(y,x) = 1 + exp(-x-y) I0(2 sqrt(xy))
    {
        const double x = 30.0, y = 22.0, s = 2.0 * std::sqrt(x * y);
        const double rhs = 1.0 + std::exp(-(std::sqrt(x) - std::sqrt(y)) * (std::sqrt(x) - std::sqrt(y))) *
                                     bed::bessel_i0_scaled(s);
        check("J(30,22) + J(22,30) identity", bed::thomas_J(x, y) + bed::thomas_J(y, x), rhs, 1e-12);
    }
    check("scaled I0 series = std::cyl_bessel_i at y = 60", bed::bessel_i0_scaled(60.0),
          std::exp(-60.0) * std::cyl_bessel_i(0.0, 60.0), 1e-15);

    std::printf("-- Linear isotherm: solver vs exact Thomas solution --\n");
    const auto L = linear_bed();
    const double tend = 3.0 * bed::tau_stoich(L);
    const double e1 = max_err_vs_thomas(L, 100, 400, tend);
    const double e2 = max_err_vs_thomas(L, 200, 800, tend);
    const double e3 = max_err_vs_thomas(L, 400, 1600, tend);
    std::printf("   max |x - J| on grids 100x400, 200x800, 400x1600: %.3e %.3e %.3e\n", e1, e2, e3);
    check("max |x - J| on the 400x1600 grid", e3, 0.0, 2e-5);
    check("observed order log2(e2/e3) (trapezoid = 2)", std::log2(e2 / e3), 2.0, 0.15);
    {
        const auto bt = bed::simulate(L, 400, tend, 1600);
        // For a linear isotherm the mean breakthrough is exactly tau_st = rho_b K L/(eps u).
        check("linear first moment = tau_st", bed::first_moment(bt), bed::tau_stoich(L),
              1e-4 * bed::tau_stoich(L));
    }

    std::printf("-- Langmuir isotherm: mass balance and constant pattern --\n");
    bed::Bed G{};
    G.length = 40.0;
    G.eps = 0.4;
    G.rho_b = 0.45;
    G.u = 20.0;
    G.k = 0.05;
    G.c0 = 0.02;
    G.iso = {100.0, 200.0, false};  // b c0 = 4, separation factor R = 0.2
    const double tst = bed::tau_stoich(G);
    const auto bt = bed::simulate(G, 2000, 1.6 * tst, 16000);
    check("first moment = tau_st (overall mass balance) [h]", bed::first_moment(bt), tst, 1e-4 * tst);
    // Constant-pattern asymptote is reached when the bed holds many MTZ lengths.
    for (double x : {0.01, 0.10, 0.50, 0.90, 0.99}) {
        char label[80];
        std::snprintf(label, sizeof label, "tau(x = %.2f) vs constant pattern [h]", x);
        check(label, bed::time_at(bt, x), bed::constant_pattern_tau(G, x), 0.5);
    }

    std::printf("-- Grid guard --\n");
    {
        // Trapezoid in z is A- but not L-stable: a cell holding more than about two
        // clean-bed transfer units lets the front ring through. simulate() refuses it.
        bool threw = false;
        try {
            bed::simulate(G, 200, tst, 1000);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        check("coarse z grid (11 transfer units per cell) rejected", threw ? 1.0 : 0.0, 1.0, 0.0);
    }

    std::printf("%s (%d failure%s)\n", failures ? "FAILED" : "ALL PASSED", failures, failures == 1 ? "" : "s");
    return failures ? 1 : 0;
}
