// Astron Built this AP
// Project 053 | Kalman filtering and sensor fusion | Li-ion state-of-charge EKF for a LEO battery
// Author: Aditya (aditya-astron)
// Date: 2026-11-21 IST
//
// Report: ZOH check, OCV observability, and the LEO Monte Carlo comparison of
// coulomb counting against three EKF variants for NMC and LFP cells.

#include <cstdio>

#include "leo.hpp"

using namespace soc;

static double rk4_zoh_error() {
    // Integrate both RC branches with RK4 (100 substeps per sample, current held)
    // over 3 orbits of the LEO profile and compare against the exact ZOH update.
    Cell c;
    Disc d = discretize(c, 1.0);
    TruthState s{0.9, 0.0, 0.0};
    double v1 = 0.0, v2 = 0.0, worst = 0.0;
    const int sub = 100;
    const double h = 1.0 / sub;
    for (int k = 0; k < static_cast<int>(3 * kOrbit); ++k) {
        double I = leo_current(k, s.z);
        auto f = [&](double v, double R, double C) { return -v / (R * C) + I / C; };
        for (int j = 0; j < sub; ++j) {
            double k1 = f(v1, c.R1, c.C1), k2 = f(v1 + 0.5 * h * k1, c.R1, c.C1);
            double k3 = f(v1 + 0.5 * h * k2, c.R1, c.C1), k4 = f(v1 + h * k3, c.R1, c.C1);
            v1 += h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
            k1 = f(v2, c.R2, c.C2), k2 = f(v2 + 0.5 * h * k1, c.R2, c.C2);
            k3 = f(v2 + 0.5 * h * k2, c.R2, c.C2), k4 = f(v2 + h * k3, c.R2, c.C2);
            v2 += h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
        }
        truth_step(s, d, I);
        worst = std::max({worst, std::fabs(v1 - s.v1), std::fabs(v2 - s.v2)});
    }
    return worst;
}

int main() {
    Cell c;
    EkfConfig cfg;
    std::printf("Cell: Q = %.1f Ah, R0 = %.0f mOhm, tau1 = %.0f s, tau2 = %.0f s\n",
                c.Q_As / 3600, c.R0 * 1e3, c.R1 * c.C1, c.R2 * c.C2);
    std::printf("Sensors: sigma_V = %.1f mV, sigma_I = %.0f mA white, bias ~ N(0, 20 mA)\n\n",
                cfg.sig_v * 1e3, cfg.sig_n * 1e3);

    std::printf("1) Exact ZOH vs RK4 (100 substeps), 3 LEO orbits: max |dv| = %.2e V\n\n",
                rk4_zoh_error());

    std::printf("2) OCV slope and the SOC resolution of one voltage sample (sigma_V / slope)\n");
    std::printf("     SOC  |   NMC mV/%%   1-sample %%  |   LFP mV/%%   1-sample %%\n");
    for (double z : {0.1, 0.3, 0.5, 0.7, 0.9}) {
        double sn = Ocv{Chem::NMC}.dv(z) * 10.0, sl = Ocv{Chem::LFP}.dv(z) * 10.0;  // mV per %
        std::printf("    %4.0f%% |  %9.2f  %10.3f   |  %9.2f  %10.3f\n", z * 100, sn,
                    cfg.sig_v * 1e3 / sn, sl, cfg.sig_v * 1e3 / sl);
    }

    std::printf("\n3) LEO cycling, 15 orbits (23.75 h), 100 runs per chemistry\n");
    std::printf("   initial SOC guess error 10 %% (1 sigma); statistics over the last orbit\n");
    for (Chem ch : {Chem::NMC, Chem::LFP}) {
        McConfig mc;
        mc.chem = ch;
        McResult r = run_leo_mc(mc);
        std::printf("\n   %s\n", chem_name(ch));
        std::printf("   method                 RMS err %%   max err %%   mean sigma %%   in 3 sigma\n");
        for (int m = 0; m < 4; ++m) {
            if (m == 0)
                std::printf("   %-21s  %9.3f  %10.3f         -            -\n",
                            method_name(Method(m)), r.rms_err[m] * 100, r.max_err[m] * 100);
            else
                std::printf("   %-21s  %9.3f  %10.3f   %10.3f   %9.1f %%\n", method_name(Method(m)),
                            r.rms_err[m] * 100, r.max_err[m] * 100, r.mean_sigz[m] * 100,
                            r.contain3[m] * 100);
        }
        for (int m = 2; m < 4; ++m)
            std::printf("   %-21s: NEES mean %.2f / per-run median %.2f (n = 4), runs > 3n: %d,\n"
                        "   %21s  worst run RMS %.3f %%, bias RMS %.2f mA, median time to\n"
                        "   %21s  stay within 1 %% SOC: %.0f s\n",
                        method_name(Method(m)), r.nees[m], r.nees_median[m], r.nees_outliers[m], "",
                        r.worst_run_err[m] * 100, r.bias_rms[m] * 1e3, "", r.t_conv_median[m]);
    }
    return 0;
}
