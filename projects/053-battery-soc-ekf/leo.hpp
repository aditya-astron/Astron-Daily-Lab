// Astron Built this AP
// Project 053 | Kalman filtering and sensor fusion | Li-ion state-of-charge EKF for a LEO battery
// Author: Aditya (aditya-astron)
// Date: 2026-11-21 IST
//
// LEO cycling scenario and Monte Carlo harness shared by the report and tests.
// 95 min orbit, 35 min eclipse. Eclipse load is 1.4 A mean with a +/-0.4 A
// payload square wave (120 s period); in sunlight the regulator charges at
// 1.0 A and tapers linearly to zero over the last 1 % below 90.5 % SOC.

#pragma once

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include "soc.hpp"

namespace soc {

constexpr double kOrbit = 5700.0, kEclipse = 2100.0;

inline double leo_current(double t, double z_true) {
    double to = std::fmod(t, kOrbit);
    if (to < kEclipse) return 1.4 + 0.4 * (std::fmod(to, 120.0) < 60.0 ? 1.0 : -1.0);
    double taper = std::clamp((0.905 - z_true) / 0.01, 0.0, 1.0);
    return -1.0 * taper;
}

enum class Method { Coulomb, Ekf3, Ekf4, Ekf4OcvInit };

inline const char* method_name(Method m) {
    switch (m) {
        case Method::Coulomb: return "coulomb count";
        case Method::Ekf3: return "EKF, no bias state";
        case Method::Ekf4: return "EKF + bias";
        default: return "EKF + bias, OCV init";
    }
}

struct McConfig {
    Chem chem = Chem::NMC;
    int runs = 100;
    int orbits = 15;
    double z0 = 0.90;
    double sig_z0 = 0.10;  // initial SOC guess error (1 sigma)
    double sig_b = 0.020;  // current-sensor bias spread [A]
    unsigned seed = 53;
};

struct McResult {
    // Statistics over the last orbit of every run.
    double rms_err[4] = {};        // SOC error RMS
    double max_err[4] = {};        // worst |SOC error|
    double contain3[4] = {};       // fraction with |e_z| <= 3 sigma_z (EKFs)
    double mean_sigz[4] = {};      // mean filter sigma_z
    // Bias-estimating filters (index 2: guess init, 3: OCV init).
    double nees[4] = {};           // mean full-state NEES (n = 4)
    double nees_median[4] = {};    // median over runs of each run's mean NEES
    int nees_outliers[4] = {};     // runs whose mean NEES exceeds 3n = 12
    double bias_rms[4] = {};       // RMS bias error at the end
    double t_conv_median[4] = {};  // median time after which |e_z| < 1 % for good [s]
    double worst_run_err[4] = {};  // largest last-orbit RMS SOC error of any run
};

inline double nees(const Ekf& f, const double* truth) {
    // e^T P^-1 e for n <= 4 via Cholesky.
    int n = f.n;
    double L[4][4] = {}, e[4];
    for (int i = 0; i < n; ++i) e[i] = f.xpost[i] - truth[i];
    for (int i = 0; i < n; ++i)
        for (int j = 0; j <= i; ++j) {
            double s = f.Ppost[i][j];
            for (int k = 0; k < j; ++k) s -= L[i][k] * L[j][k];
            L[i][j] = (i == j) ? std::sqrt(s) : s / L[j][j];
        }
    double w[4], q = 0.0;
    for (int i = 0; i < n; ++i) {
        double s = e[i];
        for (int k = 0; k < i; ++k) s -= L[i][k] * w[k];
        w[i] = s / L[i][i];
        q += w[i] * w[i];
    }
    return q;
}

inline McResult run_leo_mc(const McConfig& mc) {
    Cell cell;
    Ocv ocv{mc.chem};
    EkfConfig base;
    Disc d = discretize(cell, base.dt);
    std::mt19937_64 rng(mc.seed);
    std::normal_distribution<double> N01(0.0, 1.0);

    const int steps = static_cast<int>(mc.orbits * kOrbit / base.dt);
    const int last0 = steps - static_cast<int>(kOrbit / base.dt);
    McResult res;
    double sum_sq[4] = {}, cnt_in[4] = {}, sum_sig[4] = {};
    double sum_nees[4] = {}, sum_b2[4] = {};
    long cnt = 0;

    std::vector<double> tconv[4], run_nees[4];
    for (int r = 0; r < mc.runs; ++r) {
        double b = mc.sig_b * N01(rng);
        double zg = std::clamp(mc.z0 + mc.sig_z0 * N01(rng), 0.05, 0.99);

        EkfConfig c3 = base;
        c3.est_bias = false;
        // The bias-blind filter still knows the sensor has an offset spec, so it
        // inflates its current noise to cover sig_b; otherwise it is overconfident.
        c3.sig_n = std::sqrt(base.sig_n * base.sig_n + mc.sig_b * mc.sig_b);
        Ekf f3(cell, ocv, c3, zg, mc.sig_z0, 0.0);
        Ekf f4(cell, ocv, base, zg, mc.sig_z0, mc.sig_b);
        Ekf* f5 = nullptr;  // OCV-initialised, built on the first sample
        double zcc = zg;

        TruthState s{mc.z0, 0.0, 0.0};
        double last_bad[4] = {}, rn[4] = {}, re[4] = {};
        long rc = 0;
        for (int k = 0; k < steps; ++k) {
            double t = k * base.dt;
            double I = leo_current(t, s.z);
            double Im = I + b + base.sig_n * N01(rng);
            double y = terminal_voltage(s, cell, ocv, I) + base.sig_v * N01(rng);
            if (k == 0) {
                SocInit in = ocv_init(ocv, cell, y, Im, base.sig_v,
                                      std::hypot(base.sig_n, mc.sig_b), mc.sig_z0);
                f5 = new Ekf(cell, ocv, base, in.z, in.sigma, mc.sig_b);
            }

            f3.step(Im, y);
            f4.step(Im, y);
            f5->step(Im, y);
            const Ekf* fs[4] = {nullptr, &f3, &f4, f5};
            double ez[4] = {zcc - s.z, f3.xpost[0] - s.z, f4.xpost[0] - s.z, f5->xpost[0] - s.z};
            for (int m = 2; m < 4; ++m)
                if (std::fabs(ez[m]) > 0.01) last_bad[m] = t + base.dt;

            if (k >= last0) {
                double tr[4] = {s.z, s.v1, s.v2, b};
                for (int m = 0; m < 4; ++m) {
                    sum_sq[m] += ez[m] * ez[m];
                    res.max_err[m] = std::max(res.max_err[m], std::fabs(ez[m]));
                    if (m > 0) {
                        double sz = std::sqrt(fs[m]->Ppost[0][0]);
                        sum_sig[m] += sz;
                        if (std::fabs(ez[m]) <= 3.0 * sz) cnt_in[m] += 1.0;
                    }
                    if (m >= 2) {
                        double q = nees(*fs[m], tr);
                        sum_nees[m] += q;
                        rn[m] += q;
                        re[m] += ez[m] * ez[m];
                    }
                }
                ++cnt;
                ++rc;
            }
            zcc -= d.dtQ * Im;  // coulomb counter on the raw sensor
            truth_step(s, d, I);
        }
        const Ekf* fb[4] = {nullptr, nullptr, &f4, f5};
        for (int m = 2; m < 4; ++m) {
            sum_b2[m] += (fb[m]->bias() - b) * (fb[m]->bias() - b);
            tconv[m].push_back(last_bad[m]);
            run_nees[m].push_back(rn[m] / rc);
            if (rn[m] / rc > 12.0) ++res.nees_outliers[m];
            res.worst_run_err[m] = std::max(res.worst_run_err[m], std::sqrt(re[m] / rc));
        }
        delete f5;
    }
    for (int m = 0; m < 4; ++m) {
        res.rms_err[m] = std::sqrt(sum_sq[m] / cnt);
        res.contain3[m] = cnt_in[m] / cnt;
        res.mean_sigz[m] = sum_sig[m] / cnt;
    }
    for (int m = 2; m < 4; ++m) {
        res.nees[m] = sum_nees[m] / cnt;
        res.bias_rms[m] = std::sqrt(sum_b2[m] / mc.runs);
        std::sort(tconv[m].begin(), tconv[m].end());
        res.t_conv_median[m] = tconv[m][tconv[m].size() / 2];
        std::sort(run_nees[m].begin(), run_nees[m].end());
        res.nees_median[m] = run_nees[m][run_nees[m].size() / 2];
    }
    return res;
}

}  // namespace soc
