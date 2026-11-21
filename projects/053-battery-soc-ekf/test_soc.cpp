// Astron Built this AP
// Project 053 | Kalman filtering and sensor fusion | Li-ion state-of-charge EKF for a LEO battery
// Author: Aditya (aditya-astron)
// Date: 2026-11-21 IST
//
// Verification. The key check is (3): with a linear OCV the filter, including
// the correlated-noise time update, must reproduce the exact Gaussian
// conditional mean and covariance of x_N given y_0..y_{N-1}, computed by brute
// force from the joint distribution.

#include <cstdio>
#include <cstdlib>
#include <vector>

#include "leo.hpp"

using namespace soc;
using Vec = std::vector<double>;
using Mat = std::vector<Vec>;

static int failures = 0;

static void check(bool ok, const char* what, double val) {
    std::printf("[%s] %-62s %.3e\n", ok ? "PASS" : "FAIL", what, val);
    if (!ok) ++failures;
}

static Mat zeros(int r, int c) { return Mat(r, Vec(c, 0.0)); }

// A diag(s) B^T
static Mat adbt(const Mat& A, const Vec& s, const Mat& B) {
    Mat C = zeros(A.size(), B.size());
    for (size_t i = 0; i < A.size(); ++i)
        for (size_t j = 0; j < B.size(); ++j) {
            double t = 0.0;
            for (size_t k = 0; k < s.size(); ++k) t += A[i][k] * s[k] * B[j][k];
            C[i][j] = t;
        }
    return C;
}

// Solve S X = B for SPD S (Cholesky), B given as rows = columns of the RHS.
static Mat chol_solve(Mat S, const Mat& Bt) {
    int n = S.size();
    for (int j = 0; j < n; ++j) {
        for (int k = 0; k < j; ++k) S[j][j] -= S[j][k] * S[j][k];
        S[j][j] = std::sqrt(S[j][j]);
        for (int i = j + 1; i < n; ++i) {
            for (int k = 0; k < j; ++k) S[i][j] -= S[i][k] * S[j][k];
            S[i][j] /= S[j][j];
        }
    }
    Mat X = Bt;
    for (auto& x : X) {
        for (int i = 0; i < n; ++i) {
            for (int k = 0; k < i; ++k) x[i] -= S[i][k] * x[k];
            x[i] /= S[i][i];
        }
        for (int i = n - 1; i >= 0; --i) {
            for (int k = i + 1; k < n; ++k) x[i] -= S[k][i] * x[k];
            x[i] /= S[i][i];
        }
    }
    return X;
}

static void test_zoh() {
    Cell c;
    double worst = 0.0;
    for (double dt : {0.1, 1.0, 10.0, 60.0}) {
        Disc d = discretize(c, dt);
        double I = 2.3, v = 0.01;
        double ex = d.a1 * v + d.b1 * I;
        int sub = 20000;
        double h = dt / sub;
        auto f = [&](double x) { return -x / (c.R1 * c.C1) + I / c.C1; };
        for (int j = 0; j < sub; ++j) {
            double k1 = f(v), k2 = f(v + 0.5 * h * k1), k3 = f(v + 0.5 * h * k2), k4 = f(v + h * k3);
            v += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
        }
        worst = std::max(worst, std::fabs(v - ex));
    }
    check(worst < 1e-13, "(1) ZOH RC update vs fine RK4, dt = 0.1..60 s  [V]", worst);
}

static void test_ocv() {
    double worst = 0.0, worst_inv = 0.0;
    bool mono = true;
    for (Chem ch : {Chem::NMC, Chem::LFP}) {
        Ocv o{ch};
        double prev = -1e9;
        for (double z = 0.02; z <= 0.98; z += 0.01) {
            double h = 1e-6;
            double fd = (o.v(z + h) - o.v(z - h)) / (2 * h);
            worst = std::max(worst, std::fabs(fd - o.dv(z)) / o.dv(z));
            mono = mono && o.v(z) > prev && o.dv(z) > 0;
            prev = o.v(z);
            worst_inv = std::max(worst_inv, std::fabs(ocv_inverse(o, o.v(z)) - z));
        }
    }
    check(worst < 1e-7, "(2) dOCV/dz vs central difference, NMC and LFP  [rel]", worst);
    check(mono, "(2) OCV strictly increasing on [0.02, 0.98]", 0.0);
    check(worst_inv < 1e-12, "(2) ocv_inverse round trip", worst_inv);
}

// Brute-force conditional of x_N on y_0..y_{N-1} for the linear-OCV model.
static void test_batch_conditioning() {
    Cell c;
    EkfConfig cfg;
    cfg.dt = 30.0;
    cfg.sig_n = 0.2;   // a noisy sensor makes the cross-covariance matter
    cfg.sig_bw = 1e-4;
    const int N = 40, nx = 4, M = nx + 3 * N;
    Disc d = discretize(c, cfg.dt);
    Ocv lin{Chem::Linear};
    const double m0[4] = {0.6, 0.0, 0.0, 0.0};
    const double p0[4] = {0.05, 2e-3, 2e-3, 0.05};

    double F[4][4] = {{1, 0, 0, d.dtQ}, {0, d.a1, 0, -d.b1}, {0, 0, d.a2, -d.b2}, {0, 0, 0, 1}};
    const double G[4] = {-d.dtQ, d.b1, d.b2, 0.0};
    const double Bn[4] = {d.dtQ, -d.b1, -d.b2, 0.0};
    const double H[4] = {0.8, -1.0, -1.0, c.R0};

    // xi = [dx0 (4), n_0..n_{N-1}, e_0..e_{N-1}, w_0..w_{N-1}]
    Vec sig2(M);
    for (int i = 0; i < nx; ++i) sig2[i] = p0[i] * p0[i];
    for (int k = 0; k < N; ++k) {
        sig2[nx + k] = cfg.sig_n * cfg.sig_n;
        sig2[nx + N + k] = cfg.sig_v * cfg.sig_v;
        sig2[nx + 2 * N + k] = cfg.sig_bw * cfg.sig_bw * cfg.dt;
    }
    std::mt19937_64 rng(11);
    std::normal_distribution<double> N01(0.0, 1.0);
    Vec xi(M);
    for (int i = 0; i < M; ++i) xi[i] = std::sqrt(sig2[i]) * N01(rng);

    Vec mx(m0, m0 + 4), Im(N), my(N), y(N);
    Mat Ax = zeros(nx, M), Ay = zeros(N, M);
    for (int i = 0; i < nx; ++i) Ax[i][i] = 1.0;
    for (int k = 0; k < N; ++k) {
        Im[k] = 1.5 * std::sin(0.3 * k) + 0.5;
        my[k] = 3.2 + 0.8 * mx[0] - mx[1] - mx[2] + c.R0 * mx[3] - c.R0 * Im[k];
        for (int j = 0; j < M; ++j) {
            double t = 0.0;
            for (int i = 0; i < nx; ++i) t += H[i] * Ax[i][j];
            Ay[k][j] = t;
        }
        Ay[k][nx + k] += c.R0;
        Ay[k][nx + N + k] += 1.0;
        Vec nm(nx);
        Mat nA = zeros(nx, M);
        for (int i = 0; i < nx; ++i) {
            nm[i] = G[i] * Im[k];
            for (int l = 0; l < nx; ++l) nm[i] += F[i][l] * mx[l];
            for (int j = 0; j < M; ++j) {
                double t = 0.0;
                for (int l = 0; l < nx; ++l) t += F[i][l] * Ax[l][j];
                nA[i][j] = t;
            }
            nA[i][nx + k] += Bn[i];
        }
        nA[3][nx + 2 * N + k] += 1.0;
        mx = nm;
        Ax = nA;
    }
    // Realised measurements and truth from the sampled xi.
    for (int k = 0; k < N; ++k) {
        y[k] = my[k];
        for (int j = 0; j < M; ++j) y[k] += Ay[k][j] * xi[j];
    }
    Mat Cyy = adbt(Ay, sig2, Ay), Cxy = adbt(Ax, sig2, Ay), Cxx = adbt(Ax, sig2, Ax);
    Mat KT = chol_solve(Cyy, Cxy);  // rows: Cyy^-1 Cyx
    Vec dy(N);
    for (int k = 0; k < N; ++k) dy[k] = y[k] - my[k];
    double mean_ref[4], cov_ref[4][4];
    for (int i = 0; i < nx; ++i) {
        mean_ref[i] = mx[i];
        for (int k = 0; k < N; ++k) mean_ref[i] += KT[i][k] * dy[k];
        for (int j = 0; j < nx; ++j) {
            cov_ref[i][j] = Cxx[i][j];
            for (int k = 0; k < N; ++k) cov_ref[i][j] -= Cxy[i][k] * KT[j][k];
        }
    }

    auto run = [&](bool corr, double& dm, double& dc) {
        EkfConfig cf = cfg;
        cf.correlated = corr;
        Ekf f(c, lin, cf, m0[0], p0[0], p0[3]);
        f.P[1][1] = p0[1] * p0[1];
        f.P[2][2] = p0[2] * p0[2];
        for (int k = 0; k < N; ++k) f.step(Im[k], y[k]);
        dm = dc = 0.0;
        for (int i = 0; i < nx; ++i) {
            dm = std::max(dm, std::fabs(f.x[i] - mean_ref[i]) / std::sqrt(cov_ref[i][i]));
            for (int j = 0; j < nx; ++j)
                dc = std::max(dc, std::fabs(f.P[i][j] - cov_ref[i][j]) /
                                      std::sqrt(cov_ref[i][i] * cov_ref[j][j]));
        }
    };
    double dm, dc, dmu, dcu;
    run(true, dm, dc);
    run(false, dmu, dcu);
    check(dm < 1e-8, "(3) correlated EKF mean vs batch conditional  [sigmas]", dm);
    check(dc < 1e-8, "(3) correlated EKF covariance vs batch  [correlation units]", dc);
    check(dcu > 1e-3, "(3) ignoring E[wv] departs from the exact covariance", dcu);
    std::printf("       (uncorrelated filter: mean off by %.2e sigma)\n", dmu);
}

static void test_coulomb_drift() {
    Cell c;
    Disc d = discretize(c, 1.0);
    double b = 0.017, z = 0.9, zt = 0.9;
    int steps = 10 * static_cast<int>(kOrbit);
    TruthState s{0.9, 0, 0};
    for (int k = 0; k < steps; ++k) {
        double I = leo_current(k, s.z);
        z -= d.dtQ * (I + b);
        truth_step(s, d, I);
        zt = s.z;
    }
    double expect = -b * steps / c.Q_As;
    double err = std::fabs((z - zt) - expect);
    check(err < 1e-12, "(4) coulomb-count drift = -b t / Q over 10 orbits", err);
}

static void test_linear_nees() {
    // Linear OCV: the filter is the exact KF, so per-instant mean NEES over
    // runs must sit inside the chi-square band (n = 4, 99 %).
    Cell c;
    EkfConfig cfg;
    Ocv lin{Chem::Linear};
    Disc d = discretize(c, cfg.dt);
    const int runs = 400;
    const int probes[5] = {10, 300, 1500, 4000, 8000};
    double sum[5] = {};
    std::mt19937_64 rng(5);
    std::normal_distribution<double> N01(0.0, 1.0);
    for (int r = 0; r < runs; ++r) {
        double z0 = 0.6 + 0.1 * N01(rng), b = 0.02 * N01(rng);
        Ekf f(c, lin, cfg, 0.6, 0.1, 0.02);
        TruthState s{z0, 0.0, 0.0};
        // The filter's prior on v1, v2 is N(0, 1e-8); draw the truth from it.
        s.v1 = 1e-4 * N01(rng);
        s.v2 = 1e-4 * N01(rng);
        int p = 0;
        for (int k = 0; k <= probes[4]; ++k) {
            double I = std::sin(2e-3 * k) + 0.3;  // slow swing keeps v2 excited
            double Im = I + b + cfg.sig_n * N01(rng);
            double y = terminal_voltage(s, c, lin, I) + cfg.sig_v * N01(rng);
            f.step(Im, y);
            if (k == probes[p]) {
                double tr[4] = {s.z, s.v1, s.v2, b};
                sum[p++] += nees(f, tr);
            }
            truth_step(s, d, I);
            b += cfg.sig_bw * std::sqrt(cfg.dt) * N01(rng);
        }
    }
    // chi-square(4 runs) / runs, 99 % two-sided via Wilson-Hilferty.
    double dof = 4.0 * runs, zq = 2.5758;
    auto wh = [&](double zz) {
        double a = 2.0 / (9.0 * dof);
        return dof * std::pow(1 - a + zz * std::sqrt(a), 3) / runs;
    };
    double lo = wh(-zq), hi = wh(zq), worst = 0.0;
    bool ok = true;
    for (int p = 0; p < 5; ++p) {
        double m = sum[p] / runs;
        ok = ok && m > lo && m < hi;
        worst = std::max(worst, std::fabs(m - 4.0));
        std::printf("       step %5d: mean NEES %.3f   (99 %% band %.3f .. %.3f)\n", probes[p], m, lo,
                    hi);
    }
    check(ok, "(5) linear-OCV Monte Carlo NEES inside chi-square band", worst);
}

static void test_leo_short() {
    McConfig mc;
    mc.runs = 30;
    mc.orbits = 6;
    mc.chem = Chem::NMC;
    McResult n = run_leo_mc(mc);
    mc.chem = Chem::LFP;
    McResult l = run_leo_mc(mc);
    check(n.rms_err[0] > 0.05, "(6) coulomb counting off by > 5 % SOC after 6 orbits", n.rms_err[0]);
    check(n.contain3[1] < 0.5, "(6) bias-blind EKF overconfident (3-sigma containment < 50 %)",
          n.contain3[1]);
    check(n.rms_err[2] < 1e-3 && n.nees[2] < 6.0, "(6) NMC EKF + bias: RMS < 0.1 %, NEES < 6",
          n.rms_err[2]);
    // A run's last-orbit NEES is strongly time-correlated (the bias error barely
    // moves), so an occasional run above 3n is chance; diverged runs sit far
    // higher and in numbers.
    check(l.nees_outliers[2] >= 5, "(6) LFP EKF + bias, guess init: >= 5 of 30 runs over 3n",
          l.nees_outliers[2]);
    check(l.nees_outliers[3] <= 2 && l.nees[3] < 6.0 && l.contain3[3] > 0.95,
          "(6) LFP OCV init: <= 2 runs over 3n, NEES < 6, 3-sigma > 95 %", l.nees[3]);
}

int main() {
    test_zoh();
    test_ocv();
    test_batch_conditioning();
    test_coulomb_drift();
    test_linear_nees();
    test_leo_short();
    std::printf("\n%s (%d failure%s)\n", failures ? "FAILED" : "ALL PASSED", failures,
                failures == 1 ? "" : "s");
    return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
