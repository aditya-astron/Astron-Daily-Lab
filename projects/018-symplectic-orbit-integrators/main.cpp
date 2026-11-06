// Astron Built this AP
// Project 018 | Numerical methods | Symplectic vs. RK4 integration of long-arc Keplerian orbits
// Author: Aditya (aditya-astron)
// Date: 2026-11-06 IST
//
// Propagates a Molniya orbit (e = 0.74) for 100 revolutions with fixed-step
// RK4, leapfrog and Yoshida-4, and compares each against the exact Kepler
// solution. The point is the shape of the error, not just its size: RK4's
// energy error grows linearly with time, the symplectic methods' stays
// bounded, and all three still accumulate a secular phase (along-track) error.
#include <cmath>
#include <cstdio>
#include <vector>

#include "integrators.hpp"

using namespace orb;

namespace {

constexpr double kMu = 398600.4418;  // Earth GM [km^3/s^2]
constexpr double kDeg = M_PI / 180.0;

struct RunStats {
    double max_rel_energy;  // max over the run of |E - E0| / |E0|
    double max_rel_h;       // max over the run of |h - h0| / |h0|
    double final_pos_err;   // |r - r_exact| at the end [km]
    long force_evals;
};

RunStats run(Method m, const State& s0, double period, int steps_per_orbit, int orbits,
             std::vector<double>* energy_per_orbit = nullptr) {
    const double h = period / steps_per_orbit;
    const double E0 = energy(s0, kMu);
    const Vec3 h0 = ang_mom(s0);
    State s = s0;
    RunStats st{0.0, 0.0, 0.0, 0};
    for (int k = 0; k < orbits; ++k) {
        double orbit_max = 0.0;
        for (int i = 0; i < steps_per_orbit; ++i) {
            s = step(m, s, h, kMu);
            const double dE = std::fabs((energy(s, kMu) - E0) / E0);
            orbit_max = std::fmax(orbit_max, dE);
            st.max_rel_h = std::fmax(st.max_rel_h, norm(ang_mom(s) - h0) / norm(h0));
        }
        st.max_rel_energy = std::fmax(st.max_rel_energy, orbit_max);
        if (energy_per_orbit) energy_per_orbit->push_back(orbit_max);
    }
    const double t_end = h * steps_per_orbit * static_cast<double>(orbits);
    st.final_pos_err = norm(s.r - kepler_propagate(s0, t_end, kMu).r);
    st.force_evals = static_cast<long>(force_evals_per_step(m)) * steps_per_orbit * orbits;
    return st;
}

}  // namespace

int main() {
    const double a = 26562.0, e = 0.74;
    const State s0 = from_elements(a, e, 63.4 * kDeg, 40.0 * kDeg, 270.0 * kDeg, 0.0, kMu);
    const double period = 2.0 * M_PI * std::sqrt(a * a * a / kMu);
    const int orbits = 100;

    std::printf("Molniya orbit: a = %.0f km, e = %.2f, rp = %.0f km, ra = %.0f km, T = %.1f s\n",
                a, e, a * (1 - e), a * (1 + e), period);
    std::printf("Start at perigee, %d revolutions (%.1f days), fixed step, reference = exact Kepler\n\n",
                orbits, orbits * period / 86400.0);

    std::printf("%-9s %6s %9s %12s %12s %14s\n", "method", "N/rev", "f-evals", "max|dE/E|",
                "max|dh/h|", "final |dr| km");
    const Method methods[] = {Method::RK4, Method::Leapfrog, Method::Yoshida4};
    for (Method m : methods) {
        for (int n : {1000, 2000, 4000, 8000}) {
            const RunStats st = run(m, s0, period, n, orbits);
            std::printf("%-9s %6d %9ld %12.3e %12.3e %14.4e\n", method_name(m), n, st.force_evals,
                        st.max_rel_energy, st.max_rel_h, st.final_pos_err);
        }
    }

    // Growth of the energy error with time at a fixed step.
    std::printf("\nEnergy error vs time, N = 2000 steps/rev (max |dE/E| within the revolution)\n");
    std::printf("%-9s %12s %12s %12s %10s\n", "method", "rev 1", "rev 10", "rev 100", "100/10");
    for (Method m : methods) {
        std::vector<double> ep;
        run(m, s0, period, 2000, orbits, &ep);
        std::printf("%-9s %12.3e %12.3e %12.3e %10.2f\n", method_name(m), ep[0], ep[9], ep[99],
                    ep[99] / ep[9]);
    }

    // Equal-cost comparison: 12000 force evaluations per revolution.
    std::printf("\nEqual cost, 12000 force evaluations per revolution, %d revs\n", orbits);
    std::printf("%-9s %6s %12s %14s\n", "method", "N/rev", "max|dE/E|", "final |dr| km");
    for (Method m : methods) {
        const int n = 12000 / force_evals_per_step(m);
        const RunStats st = run(m, s0, period, n, orbits);
        std::printf("%-9s %6d %12.3e %14.4e\n", method_name(m), n, st.max_rel_energy,
                    st.final_pos_err);
    }
    return 0;
}
