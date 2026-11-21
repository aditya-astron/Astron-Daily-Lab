// Astron Built this AP
// Project 053 | Kalman filtering and sensor fusion | Li-ion state-of-charge EKF for a LEO battery
// Author: Aditya (aditya-astron)
// Date: 2026-11-21 IST
//
// Second-order Thevenin cell model, exact zero-order-hold discretisation, and
// an EKF for state of charge with an optional current-sensor bias state.
//
// Sign convention: current I > 0 discharges the cell.
// The current sensor reads Im = I + b + n, with constant bias b and white
// noise n. The same n drives the state (through the coulomb count and the RC
// branches) and the voltage residual (through R0), so process and measurement
// noise are correlated. The filter can carry that cross-covariance exactly.

#pragma once

#include <algorithm>
#include <cmath>

namespace soc {

enum class Chem { NMC, LFP, Linear };

inline const char* chem_name(Chem c) {
    switch (c) {
        case Chem::NMC: return "NMC";
        case Chem::LFP: return "LFP";
        default: return "linear";
    }
}

// Open-circuit voltage vs. state of charge z in [0, 1].
// NMC: sloped curve, ~11 mV per % SOC mid-range, log terms give the knees.
// LFP: flat two-phase plateau (~0.5 mV per % SOC) with steep ends.
// Linear: used for exact linear-Gaussian checks; never clamped.
struct Ocv {
    Chem chem;
    static constexpr double zmin = 0.005, zmax = 0.995;

    double v(double z) const {
        if (chem == Chem::Linear) return 3.2 + 0.8 * z;
        z = std::clamp(z, zmin, zmax);
        if (chem == Chem::NMC)
            return 3.0 + 1.2 * z - 0.9 * z * z + 0.9 * z * z * z + 0.04 * std::log(z) -
                   0.03 * std::log(1.0 - z);
        return 3.28 + 0.05 * (z - 0.5) - 0.35 * std::exp(-z / 0.06) +
               0.20 * std::exp(-(1.0 - z) / 0.04);
    }

    double dv(double z) const {
        if (chem == Chem::Linear) return 0.8;
        if (z < zmin || z > zmax) return 0.0;
        if (chem == Chem::NMC)
            return 1.2 - 1.8 * z + 2.7 * z * z + 0.04 / z + 0.03 / (1.0 - z);
        return 0.05 + (0.35 / 0.06) * std::exp(-z / 0.06) +
               (0.20 / 0.04) * std::exp(-(1.0 - z) / 0.04);
    }
};

// One 3.0 Ah 18650-class cell. tau1 = 30 s (charge transfer), tau2 = 1000 s
// (diffusion).
struct Cell {
    double Q_As = 3.0 * 3600.0;
    double R0 = 0.030;
    double R1 = 0.015, C1 = 2000.0;
    double R2 = 0.020, C2 = 50000.0;
};

// Exact ZOH: with I constant over dt, dv/dt = -v/(RC) + I/C integrates to
// v+ = a v + R (1 - a) I, a = exp(-dt / RC). No integration error at any dt.
struct Disc {
    double a1, a2, b1, b2, dtQ;
};

inline Disc discretize(const Cell& c, double dt) {
    Disc d;
    d.a1 = std::exp(-dt / (c.R1 * c.C1));
    d.a2 = std::exp(-dt / (c.R2 * c.C2));
    d.b1 = c.R1 * (1.0 - d.a1);
    d.b2 = c.R2 * (1.0 - d.a2);
    d.dtQ = dt / c.Q_As;
    return d;
}

struct TruthState {
    double z, v1, v2;
};

inline void truth_step(TruthState& s, const Disc& d, double I) {
    s.z -= d.dtQ * I;
    s.v1 = d.a1 * s.v1 + d.b1 * I;
    s.v2 = d.a2 * s.v2 + d.b2 * I;
}

inline double terminal_voltage(const TruthState& s, const Cell& c, const Ocv& ocv, double I) {
    return ocv.v(s.z) - s.v1 - s.v2 - c.R0 * I;
}

// Invert the OCV curve by bisection (it is strictly increasing).
inline double ocv_inverse(const Ocv& ocv, double v) {
    double lo = Ocv::zmin, hi = Ocv::zmax;
    for (int i = 0; i < 60; ++i) {
        double mid = 0.5 * (lo + hi);
        (ocv.v(mid) < v ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
}

struct SocInit {
    double z, sigma;
};

// First-sample initialisation with the RC branches relaxed: OCV = y + R0 Im.
// The sigma maps voltage noise and the current-sensor spread through R0 and
// the local slope, capped at the a-priori spread.
inline SocInit ocv_init(const Ocv& ocv, const Cell& c, double y0, double Im0, double sig_v,
                        double sig_i, double sig_cap) {
    double z = ocv_inverse(ocv, y0 + c.R0 * Im0);
    double sv = std::sqrt(sig_v * sig_v + c.R0 * c.R0 * sig_i * sig_i);
    double slope = std::max(ocv.dv(z), 1e-9);
    return {z, std::min(sig_cap, sv / slope)};
}

struct EkfConfig {
    bool est_bias = true;     // augment with current-sensor bias b
    bool correlated = true;   // use E[w v] = Bn R0 sig_n^2 in the time update
    double dt = 1.0;
    double sig_v = 2e-3;      // voltage ADC noise [V]
    double sig_n = 10e-3;     // current sensor white noise [A]
    double sig_bw = 1e-6;     // bias random-walk density [A/sqrt(s)]
};

// State x = [z, v1, v2, b]; n = 3 drops b.
class Ekf {
public:
    int n;
    double x[4] = {0, 0, 0, 0};
    double P[4][4] = {};
    double xpost[4] = {0, 0, 0, 0};  // x_{k|k} from the last step
    double Ppost[4][4] = {};
    double nis = 0.0;

    Ekf(const Cell& cell, Ocv ocv, EkfConfig cfg, double z0, double sig_z0, double sig_b0)
        : n(cfg.est_bias ? 4 : 3), cell_(cell), ocv_(ocv), cfg_(cfg), d_(discretize(cell, cfg.dt)) {
        x[0] = z0;
        P[0][0] = sig_z0 * sig_z0;
        P[1][1] = P[2][2] = 1e-8;  // RC branches start relaxed, known to ~0.1 mV
        if (n == 4) P[3][3] = sig_b0 * sig_b0;

        double bn[4] = {d_.dtQ, -d_.b1, -d_.b2, 0.0};
        double s2 = cfg_.sig_n * cfg_.sig_n;
        R_ = cell_.R0 * cell_.R0 * s2 + cfg_.sig_v * cfg_.sig_v;
        for (int i = 0; i < 4; ++i) {
            S_[i] = bn[i] * cell_.R0 * s2;
            for (int j = 0; j < 4; ++j) Q_[i][j] = bn[i] * bn[j] * s2;
        }
        Q_[3][3] += cfg_.sig_bw * cfg_.sig_bw * cfg_.dt;

        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) F_[i][j] = 0.0;
        F_[0][0] = 1.0;
        F_[1][1] = d_.a1;
        F_[2][2] = d_.a2;
        F_[3][3] = 1.0;
        if (n == 4) {  // states see the corrected current Im - b
            F_[0][3] = d_.dtQ;
            F_[1][3] = -d_.b1;
            F_[2][3] = -d_.b2;
        }
    }

    double bias() const { return n == 4 ? x[3] : 0.0; }

    double predict_y(const double* s, double Im) const {
        double b = n == 4 ? s[3] : 0.0;
        return ocv_.v(s[0]) - s[1] - s[2] - cell_.R0 * (Im - b);
    }

    // One sample: measurement update with y_k, then time update to k+1.
    void step(double Im, double y) {
        double H[4];
        jacobian(x, H);
        double PH[4];
        for (int i = 0; i < n; ++i) {
            PH[i] = 0.0;
            for (int j = 0; j < n; ++j) PH[i] += P[i][j] * H[j];
        }
        double Sv = R_;
        for (int i = 0; i < n; ++i) Sv += H[i] * PH[i];
        double nu = y - predict_y(x, Im);
        nis = nu * nu / Sv;

        double K[4];
        for (int i = 0; i < n; ++i) {
            K[i] = PH[i] / Sv;
            x[i] += K[i] * nu;
        }
        // Joseph form keeps P symmetric positive definite.
        double A[4][4], T[4][4];
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) A[i][j] = (i == j ? 1.0 : 0.0) - K[i] * H[j];
        mul(A, P, T);
        mul_bt(T, A, P);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) P[i][j] += K[i] * K[j] * R_;

        for (int i = 0; i < n; ++i) {
            xpost[i] = x[i];
            for (int j = 0; j < n; ++j) Ppost[i][j] = P[i][j];
        }

        // Time update. With E[w_k v_k] = S, write w_k = S v_k / R + w~_k where
        // v_k = y_k - h(x_k); then x_{k+1} = f(x_k) + S (y_k - h(x_k)) / R + w~_k,
        // so F~ = F - S H / R and Q~ = Q - S S^T / R.
        double b = n == 4 ? x[3] : 0.0;
        double Icorr = Im - b;
        double gain[4] = {0, 0, 0, 0};
        double Ft[4][4];
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) Ft[i][j] = F_[i][j];
        double Qt[4][4];
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) Qt[i][j] = Q_[i][j];
        if (cfg_.correlated) {
            double r = y - predict_y(x, Im);
            jacobian(x, H);
            for (int i = 0; i < n; ++i) {
                gain[i] = S_[i] * r / R_;
                for (int j = 0; j < n; ++j) {
                    Ft[i][j] -= S_[i] * H[j] / R_;
                    Qt[i][j] -= S_[i] * S_[j] / R_;
                }
            }
        }
        x[0] = x[0] - d_.dtQ * Icorr + gain[0];
        x[1] = d_.a1 * x[1] + d_.b1 * Icorr + gain[1];
        x[2] = d_.a2 * x[2] + d_.b2 * Icorr + gain[2];
        if (n == 4) x[3] += gain[3];

        mul(Ft, P, T);
        mul_bt(T, Ft, P);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) P[i][j] += Qt[i][j];
    }

private:
    Cell cell_;
    Ocv ocv_;
    EkfConfig cfg_;
    Disc d_;
    double R_;
    double S_[4];
    double Q_[4][4];
    double F_[4][4];

    void jacobian(const double* s, double* H) const {
        H[0] = ocv_.dv(s[0]);
        H[1] = -1.0;
        H[2] = -1.0;
        H[3] = cell_.R0;
    }
    void mul(const double A[4][4], const double B[4][4], double C[4][4]) const {
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                double s = 0.0;
                for (int k = 0; k < n; ++k) s += A[i][k] * B[k][j];
                C[i][j] = s;
            }
    }
    void mul_bt(const double A[4][4], const double B[4][4], double C[4][4]) const {
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                double s = 0.0;
                for (int k = 0; k < n; ++k) s += A[i][k] * B[j][k];
                C[i][j] = s;
            }
    }
};

}  // namespace soc
