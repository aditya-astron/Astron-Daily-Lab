// Astron Built this AP
// Project 018 | Numerical methods | Symplectic vs. RK4 integration of long-arc Keplerian orbits
// Author: Aditya (aditya-astron)
// Date: 2026-11-06 IST
#include <cmath>
#include <cstdio>
#include <vector>

#include "integrators.hpp"

using namespace orb;

namespace {

constexpr double kMu = 398600.4418;
constexpr double kDeg = M_PI / 180.0;
int g_fail = 0;

void check(bool ok, const char* what, double value) {
    std::printf("[%s] %-62s %.6g\n", ok ? "PASS" : "FAIL", what, value);
    if (!ok) ++g_fail;
}

double period_of(double a) { return 2.0 * M_PI * std::sqrt(a * a * a / kMu); }

State propagate(Method m, State s, double h, long n) {
    for (long i = 0; i < n; ++i) s = step(m, s, h, kMu);
    return s;
}

}  // namespace

int main() {
    // --- Reference solution -------------------------------------------------
    // Hand check of from_elements at perigee: r = a(1-e), v = sqrt(mu (1+e)/(a(1-e))).
    {
        const double a = 26562.0, e = 0.74;
        const State s = from_elements(a, e, 63.4 * kDeg, 40 * kDeg, 270 * kDeg, 0.0, kMu);
        const double vp = std::sqrt(kMu * (1 + e) / (a * (1 - e)));
        check(std::fabs(norm(s.r) - a * (1 - e)) < 1e-9, "perigee radius a(1-e) [km err]",
              norm(s.r) - a * (1 - e));
        check(std::fabs(norm(s.v) - vp) < 1e-12, "perigee speed from vis-viva [km/s err]",
              norm(s.v) - vp);
        check(std::fabs(dot(s.r, s.v)) < 1e-6, "r.v = 0 at perigee", dot(s.r, s.v));
        // Molniya: inclination recovered from h_z / |h|
        const Vec3 hv = ang_mom(s);
        check(std::fabs(std::acos(hv[2] / norm(hv)) / kDeg - 63.4) < 1e-10,
              "inclination from angular momentum [deg]", std::acos(hv[2] / norm(hv)) / kDeg);
    }
    {
        const double a = 26562.0, e = 0.74;
        const State s0 = from_elements(a, e, 63.4 * kDeg, 40 * kDeg, 270 * kDeg, 1.1, kMu);
        const double T = period_of(a);
        const State sT = kepler_propagate(s0, T, kMu);
        check(norm(sT.r - s0.r) < 1e-6, "Kepler: full period returns to start [km]", norm(sT.r - s0.r));
        const State s37 = kepler_propagate(s0, 37.3 * T, kMu);
        const State sback = kepler_propagate(s37, -37.3 * T, kMu);
        check(norm(sback.r - s0.r) < 1e-5, "Kepler: 37.3 revs forward then back [km]",
              norm(sback.r - s0.r));
        // Same point reached by mean anomaly directly.
        const double n = std::sqrt(kMu / (a * a * a));
        const State sm = from_elements(a, e, 63.4 * kDeg, 40 * kDeg, 270 * kDeg, 1.1 + n * 5000.0, kMu);
        const State sk = kepler_propagate(s0, 5000.0, kMu);
        check(norm(sk.r - sm.r) < 1e-7, "Kepler f,g vs. elements at M0 + n dt [km]", norm(sk.r - sm.r));
        check(std::fabs((energy(s37, kMu) - energy(s0, kMu)) / energy(s0, kMu)) < 1e-13,
              "Kepler: energy invariant", energy(s37, kMu) - energy(s0, kMu));
    }

    // --- Order of accuracy --------------------------------------------------
    // One revolution of an e = 0.3 orbit at N and 2N steps; the error ratio
    // must approach 2^p.
    {
        const double a = 9000.0, e = 0.3;
        const State s0 = from_elements(a, e, 30 * kDeg, 10 * kDeg, 50 * kDeg, 0.0, kMu);
        const double T = period_of(a);
        const State exact = kepler_propagate(s0, T, kMu);
        struct Expect { Method m; double p; long n; };
        for (Expect ex : {Expect{Method::RK4, 4.0, 1600}, Expect{Method::Leapfrog, 2.0, 2000},
                          Expect{Method::Yoshida4, 4.0, 400}}) {
            const double e1 = norm(propagate(ex.m, s0, T / ex.n, ex.n).r - exact.r);
            const double e2 = norm(propagate(ex.m, s0, T / (2 * ex.n), 2 * ex.n).r - exact.r);
            const double p = std::log2(e1 / e2);
            char label[96];
            std::snprintf(label, sizeof label, "%s observed order (expect %.0f)", method_name(ex.m), ex.p);
            check(std::fabs(p - ex.p) < 0.1, label, p);
        }
        // Yoshida coefficients: consistency w1 + w0 + w1 = 1 and the
        // third-order condition 2 w1^3 + w0^3 = 0.
        const double c = std::cbrt(2.0), w1 = 1 / (2 - c), w0 = -c / (2 - c);
        check(std::fabs(2 * w1 + w0 - 1) < 1e-15, "Yoshida sum of weights - 1", 2 * w1 + w0 - 1);
        check(std::fabs(2 * w1 * w1 * w1 + w0 * w0 * w0) < 1e-14, "Yoshida 2 w1^3 + w0^3", 2 * w1 * w1 * w1 + w0 * w0 * w0);
    }

    // --- Geometric properties -----------------------------------------------
    const double a = 26562.0, e = 0.74;
    const State s0 = from_elements(a, e, 63.4 * kDeg, 40 * kDeg, 270 * kDeg, 0.0, kMu);
    const double T = period_of(a);
    const double E0 = energy(s0, kMu);
    const double h0 = norm(ang_mom(s0));
    const int N = 2000, revs = 100;
    const double h = T / N;

    // Per revolution: max |dE/E| within the revolution, and |dE/E| sampled
    // at the end of each revolution (same orbital phase, so the periodic
    // perigee excursion drops out and only the secular part remains).
    std::vector<double> rk_end;
    auto envelope = [&](Method m, std::vector<double>* at_end) {
        std::vector<double> env;
        State s = s0;
        double hdrift = 0.0;
        for (int k = 0; k < revs; ++k) {
            double mx = 0.0;
            for (int i = 0; i < N; ++i) {
                s = step(m, s, h, kMu);
                mx = std::fmax(mx, std::fabs((energy(s, kMu) - E0) / E0));
                hdrift = std::fmax(hdrift, norm(ang_mom(s) - ang_mom(s0)) / h0);
            }
            env.push_back(mx);
            if (at_end) at_end->push_back(std::fabs((energy(s, kMu) - E0) / E0));
        }
        env.push_back(hdrift);  // last entry: max angular momentum drift
        return env;
    };

    const std::vector<double> rk = envelope(Method::RK4, &rk_end);
    const std::vector<double> lf = envelope(Method::Leapfrog, nullptr);
    const std::vector<double> y4 = envelope(Method::Yoshida4, nullptr);

    // RK4: secular energy error grows linearly, so end-of-rev 100 / rev 10 = 10.
    check(std::fabs(rk_end[99] / rk_end[9] - 10.0) < 0.1, "RK4 end-of-rev energy error rev100/rev10 (linear: 10)",
          rk_end[99] / rk_end[9]);
    // Symplectic: bounded, same envelope at rev 1 and rev 100.
    check(lf[99] / lf[0] < 1.05, "Leapfrog energy envelope rev100/rev1 (bounded)", lf[99] / lf[0]);
    check(y4[99] / y4[0] < 1.05, "Yoshida4 energy envelope rev100/rev1 (bounded)", y4[99] / y4[0]);
    // Angular momentum: exact to round-off for the symplectic splittings.
    check(lf[revs] < 1e-12, "Leapfrog |dh|/h over 100 revs (round-off only)", lf[revs]);
    check(y4[revs] < 1e-12, "Yoshida4 |dh|/h over 100 revs (round-off only)", y4[revs]);
    check(rk[revs] > 1e-9, "RK4 |dh|/h over 100 revs (truncation drift)", rk[revs]);

    // Time reversibility: leapfrog forward, flip velocity, forward again.
    {
        const State f = propagate(Method::Leapfrog, s0, h, 10L * N);
        const State back = propagate(Method::Leapfrog, {f.r, -1.0 * f.v}, h, 10L * N);
        check(norm(back.r - s0.r) < 1e-6, "Leapfrog reversibility over 10 revs [km]", norm(back.r - s0.r));
        const State fr = propagate(Method::RK4, s0, h, 10L * N);
        const State br = propagate(Method::RK4, {fr.r, -1.0 * fr.v}, h, 10L * N);
        check(norm(br.r - s0.r) > 1e3 * norm(back.r - s0.r), "RK4 is not reversible (error >> leapfrog)",
              norm(br.r - s0.r));
    }

    std::printf("\n%s (%d failure%s)\n", g_fail ? "FAILED" : "ALL PASSED", g_fail, g_fail == 1 ? "" : "s");
    return g_fail ? 1 : 0;
}
