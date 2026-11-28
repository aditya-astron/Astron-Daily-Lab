// Astron Built this AP
// Project 068 | Numerical methods / Orbital mechanics | Spherical-harmonic gravity and GEO longitude drift
// Author: Aditya (aditya-astron)
// Date: 2026-11-28 IST
//
// Study: where does the tesseral field park an uncontrolled geostationary
// satellite, how fast does it librate, and what does east-west station
// keeping cost at a given slot?

#include <cstdio>

#include "geo.hpp"

using namespace geo;

int main() {
    const Field egm = grav::egm96_4x4();
    const double r = sync_radius(egm);

    std::printf("EGM96 4x4, GM = %.9e m^3/s^2, R = %.1f m, synchronous radius %.3f km\n\n",
                egm.GM, egm.R, r / 1e3);
    std::printf("J2 = %.6e  J3 = %.6e  J4 = %.6e  lambda22 = %.3f deg\n\n",
                -std::sqrt(5.0) * egm.c(2, 0), -std::sqrt(7.0) * egm.c(3, 0), -3.0 * egm.c(4, 0),
                0.5 * std::atan2(egm.s(2, 2), egm.c(2, 2)) / kDeg);

    // 1. Equilibria, degree by degree
    for (int n : {2, 3, 4}) {
        Field f = egm.truncated(n);
        std::printf("Equilibria, field to degree %d:\n", n);
        for (auto& e : equilibria(f, r))
            std::printf("  %8.3f E  (%7.3f %s)  %-8s  d(lam_ddot)/d(lam) = %+.3e 1/day^2\n",
                        e.lam_deg, e.lam_deg > 180 ? 360 - e.lam_deg : e.lam_deg,
                        e.lam_deg > 180 ? "W" : "E", e.stable ? "stable" : "unstable",
                        e.curvature * kDay * kDay);
        std::printf("\n");
    }

    // 2. Drift acceleration and station-keeping cost around the equator
    std::printf("  lon (E)   lam_ddot (deg/day^2)   a_T (m/s^2)   EW dv (m/s/yr)\n");
    double peak = 0.0, peak_lon = 0.0;
    for (int d = 0; d < 360; ++d) {
        double a = along_track(egm, d * kDeg, r);
        if (std::fabs(a) > peak) { peak = std::fabs(a); peak_lon = d; }
        if (d % 30 == 0)
            std::printf("  %6d    %+12.6f        %+.3e     %6.3f\n", d,
                        lam_ddot(egm, d * kDeg, r) / kDeg * kDay * kDay, a, std::fabs(a) * kYear);
    }
    std::printf("  peak |a_T| %.3e m/s^2 at %.0f E -> %.3f m/s/yr\n\n", peak, peak_lon, peak * kYear);

    // 3. Free libration about the stable points: 1D reduction vs full 3D orbit
    std::printf("Free libration, released at rest (period in days):\n");
    std::printf("  well      release      1D model    3D orbit    diff\n");
    for (auto& e : equilibria(egm, r)) {
        if (!e.stable) continue;
        double T0 = 2 * kPi / std::sqrt(-e.curvature) / kDay;
        std::printf("  %7.2f E  small-angle  %9.1f\n", e.lam_deg, T0);
        for (double off : {1.0, 20.0, 60.0}) {
            double l0 = (e.lam_deg + off) * kDeg;
            double T1 = libration_period_1d(egm, r, l0, e.lam_deg * kDeg) / kDay;
            Libration3D L = libration_period_3d(egm, l0, e.lam_deg * kDeg, 2);
            std::printf("            +%4.0f deg    %9.1f   %9.1f   %+.3f %%   (swing %+.1f/%+.1f deg)\n",
                        off, T1, L.period / kDay, 100 * (L.period / kDay - T1) / T1, L.min_deg, L.max_deg);
        }
    }
    std::printf("\n");

    // 4. East-west station keeping in a +-0.05 deg box, 3D simulation
    std::printf("East-west box keeping, +-0.05 deg, paired tangential burns, 3 years:\n");
    std::printf("  slot    predicted |a_T|*yr   simulated     cycle     burns   max excursion\n");
    for (double slot : {30.0, 60.0, 100.0, 200.0}) {
        double pred = std::fabs(along_track(egm, slot * kDeg, r)) * kYear;
        StationKeeping sk = east_west_keeping(egm, slot * kDeg, 0.05 * kDeg, 3.0);
        std::printf("  %5.0f E    %7.3f m/s/yr     %7.3f    %6.1f d    %4d     %.4f deg\n",
                    slot, pred, sk.dv_rate, sk.cycle_days, sk.burns, sk.max_excursion);
    }
    return 0;
}
