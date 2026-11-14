// Astron Built this AP
// Project 037 | ECLSS | Water recovery sorbent bed breakthrough and LUB scale-up
// Author: Aditya (aditya-astron)
// Date: 2026-11-14 IST
//
// Sizing the activated-carbon polishing bed of a 4-crew water processor.
// A 20 cm pilot bed is simulated, its length of unused bed (LUB) is used to
// scale the bed to a one-year life at a 1 % breakthrough limit, and the scaled
// bed is simulated to see how well the shortcut held.

#include <cmath>
#include <cstdio>

#include "sorbent_bed.hpp"

namespace {

constexpr double kFlow = 600.0;        // processed water [cm^3/h] (14.4 kg/day)
constexpr double kDiameter = 10.0;     // bed inner diameter [cm]
constexpr double kLimit = 0.01;        // allowed outlet fraction c/c0
constexpr double kLifeTarget = 8760.0; // required life [h] (one year)
constexpr double kDz = 0.01;           // z step [cm]: 0.06 clean-bed transfer units per cell
constexpr double kDtau = 1.0;          // tau step [h]

bed::Bed make_bed(double length) {
    bed::Bed B{};
    const double area = M_PI * kDiameter * kDiameter / 4.0;
    B.length = length;
    B.eps = 0.40;
    B.rho_b = 0.45;                    // g/cm^3, granular activated carbon
    B.u = kFlow / (area * B.eps);      // interstitial velocity
    B.k = 0.005;                       // LDF constant, intraparticle-limited [1/h]
    B.c0 = 20e-3;                      // 20 mg/L model organic in the feed
    B.iso = {100.0, 200.0, false};     // qm = 100 mg/g, b = 0.2 L/mg  -> R = 0.2
    return B;
}

bed::Breakthrough run(const bed::Bed& B, double tau_end) {
    const int nz = static_cast<int>(std::ceil(B.length / kDz));
    const int nt = static_cast<int>(std::ceil(tau_end / kDtau));
    return bed::simulate(B, nz, tau_end, nt);
}

double sorbent_kg(const bed::Bed& B) {
    return B.rho_b * M_PI * kDiameter * kDiameter / 4.0 * B.length / 1000.0;
}

}  // namespace

int main() {
    const bed::Bed P = make_bed(20.0);
    const double R = 1.0 / (1.0 + P.iso.b * P.c0);
    const double q0 = P.iso.q(P.c0);
    const double transit = P.length / P.u;

    std::printf("Water processor polishing bed, model organic at %.0f mg/L\n", P.c0 * 1e3);
    std::printf("  flow %.0f mL/h, ID %.0f cm, u = %.2f cm/h, q0 = %.1f mg/g, R = %.2f\n", kFlow, kDiameter,
                P.u, q0, R);
    std::printf("  wave speed %.4f cm/day, clean-bed rate %.2f transfer units/cm\n\n",
                24.0 * P.length / bed::tau_stoich(P),
                P.rho_b * P.k * P.iso.dq(0.0) / (P.eps * P.u));

    // Pilot bed
    const double tstP = bed::tau_stoich(P);
    const auto btP = run(P, 1.5 * tstP);
    const double tbP = bed::time_at(btP, kLimit);
    const double lubP = P.length * (1.0 - tbP / tstP);
    std::printf("Pilot bed L = %.1f cm (%.2f kg carbon)\n", P.length, sorbent_kg(P));
    std::printf("  stoichiometric time   %8.1f h  (%.1f days)\n", tstP + transit, (tstP + transit) / 24.0);
    std::printf("  1%% breakthrough       %8.1f h  (%.1f days)\n", tbP + transit, (tbP + transit) / 24.0);
    std::printf("  50%% / 99%%             %8.1f h / %.1f h\n", bed::time_at(btP, 0.5) + transit,
                bed::time_at(btP, 0.99) + transit);
    std::printf("  utilisation at limit  %8.3f\n", tbP / tstP);
    std::printf("  LUB                   %8.3f cm\n", lubP);

    // Asymptotic (constant-pattern) LUB for comparison: depends only on k, R and the wave speed.
    const double lubCP = (tstP - bed::constant_pattern_tau(P, kLimit)) * P.length / tstP;
    std::printf("  constant-pattern LUB  %8.3f cm\n\n", lubCP);

    std::printf("  outlet history (days : c/c0)\n");
    for (double d = 160.0; d <= 260.0; d += 10.0) {
        const size_t j = static_cast<size_t>(std::lround((d * 24.0 - transit) / kDtau));
        if (j < btP.x.size()) std::printf("    %5.0f : %.4f\n", d, btP.x[j]);
    }

    // LUB scale-up: equilibrium length for the target life plus the unused zone.
    const double wave = P.length / tstP;  // cm per hour of retarded time
    const double Leq = wave * kLifeTarget;
    const double Lsized = Leq + lubP;
    std::printf("\nSizing for %.0f days at c/c0 <= %.2f\n", kLifeTarget / 24.0, kLimit);
    std::printf("  equilibrium length (no MTZ)  %7.2f cm\n", Leq);
    std::printf("  LUB scale-up length          %7.2f cm\n", Lsized);

    const bed::Bed S = make_bed(Lsized);
    const auto btS = run(S, 1.15 * bed::tau_stoich(S));
    const double tbS = bed::time_at(btS, kLimit) + S.length / S.u;
    std::printf("  simulated life of that bed   %7.1f days (target %.0f, error %+.2f %%)\n", tbS / 24.0,
                kLifeTarget / 24.0, 100.0 * (tbS - kLifeTarget) / kLifeTarget);
    std::printf("  carbon mass                  %7.2f kg\n", sorbent_kg(S));

    const bed::Bed E = make_bed(Leq);
    const auto btE = run(E, 1.15 * bed::tau_stoich(E));
    const double tbE = bed::time_at(btE, kLimit) + E.length / E.u;
    std::printf("  equilibrium-sized bed lasts  %7.1f days (%.1f days short)\n", tbE / 24.0,
                (kLifeTarget - tbE) / 24.0);
    return 0;
}
