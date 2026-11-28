// Astron Built this AP
// Project 068 | Numerical methods / Orbital mechanics | Spherical-harmonic gravity and GEO longitude drift
// Author: Aditya (aditya-astron)
// Date: 2026-11-28 IST
//
// Verification: Legendre closed forms and orthonormality, two independent
// acceleration paths against each other and against a numerical gradient of
// the potential, J2 closed form, degree-2 equilibria by hand, published
// equilibrium longitudes, 1D drift model against the full 3D orbit, and the
// station-keeping cost law.

#include <cstdio>
#include <random>

#include "geo.hpp"

using namespace geo;

static int failures = 0;

static void check(bool ok, const char* what, double got, double want) {
    std::printf("  [%s] %-58s got %+.6e  ref %+.6e\n", ok ? "pass" : "FAIL", what, got, want);
    if (!ok) ++failures;
}

static void check_rel(const char* what, double got, double want, double tol) {
    check(std::fabs(got - want) <= tol * std::fabs(want), what, got, want);
}

static double vnorm(const Vec3& a) { return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]); }
static double vdiff(const Vec3& a, const Vec3& b) {
    return vnorm({a[0] - b[0], a[1] - b[1], a[2] - b[2]});
}

// Gauss-Legendre nodes and weights on [-1, 1] by Newton on P_n.
static void gauss_legendre(int n, std::vector<double>& x, std::vector<double>& w) {
    x.resize(n); w.resize(n);
    for (int i = 0; i < n; ++i) {
        double z = std::cos(kPi * (i + 0.75) / (n + 0.5)), dp = 0.0;
        for (int it = 0; it < 100; ++it) {
            double p0 = 1.0, p1 = z;
            for (int k = 2; k <= n; ++k) { double p2 = ((2 * k - 1) * z * p1 - (k - 1) * p0) / k; p0 = p1; p1 = p2; }
            dp = n * (z * p1 - p0) / (z * z - 1);
            double dz = p1 / dp;
            z -= dz;
            if (std::fabs(dz) < 1e-16) break;
        }
        x[i] = z;
        w[i] = 2.0 / ((1 - z * z) * dp * dp);
    }
}

static Vec3 fd_gradient(const Field& f, const Vec3& x, double h) {
    Vec3 g{};
    for (int i = 0; i < 3; ++i) {
        Vec3 p = x, m = x;
        p[i] += h; m[i] -= h;
        Vec3 p2 = x, m2 = x;
        p2[i] += 2 * h; m2[i] -= 2 * h;
        // fourth-order central difference
        g[i] = (8 * (grav::potential(f, p) - grav::potential(f, m))
                - (grav::potential(f, p2) - grav::potential(f, m2))) / (12 * h);
    }
    return g;
}

int main() {
    const Field egm = grav::egm96_4x4();
    const double r = sync_radius(egm);

    std::printf("Coefficients\n");
    check_rel("J2 = -sqrt(5) C20", -std::sqrt(5.0) * egm.c(2, 0), 1.08263e-3, 1e-5);
    check_rel("J3 = -sqrt(7) C30", -std::sqrt(7.0) * egm.c(3, 0), -2.5327e-6, 1e-4);
    check_rel("J4 = -3 C40", -3.0 * egm.c(4, 0), -1.6196e-6, 1e-4);
    check_rel("lambda22 = atan2(S22, C22)/2 (deg)", 0.5 * std::atan2(egm.s(2, 2), egm.c(2, 2)) / kDeg, -14.93, 1e-3);
    check_rel("synchronous radius (km)", r / 1e3, 42164.17, 1e-6);

    std::printf("Normalized Legendre functions\n");
    {
        double worst = 0.0;
        for (double phi : {-1.2, -0.4, 0.0, 0.3, 1.1}) {
            double t = std::sin(phi), u = std::cos(phi);
            std::vector<double> P;
            grav::legendre(4, t, u, P);
            auto at = [&](int n, int m) { return P[Field::tri(n, m)]; };
            double ref[][3] = {
                {2, 0, std::sqrt(5.0) * (3 * t * t - 1) / 2},
                {2, 1, std::sqrt(15.0) * t * u},
                {2, 2, std::sqrt(15.0) / 2 * u * u},
                {3, 0, std::sqrt(7.0) * (5 * t * t * t - 3 * t) / 2},
                {3, 3, std::sqrt(35.0 / 8) * u * u * u},
                {4, 0, 3.0 * (35 * t * t * t * t - 30 * t * t + 3) / 8},
            };
            for (auto& rw : ref) worst = std::max(worst, std::fabs(at(int(rw[0]), int(rw[1])) - rw[2]));
        }
        check(worst < 1e-14, "closed forms Pbar20..Pbar40, max abs error", worst, 0.0);

        // (1/4pi) integral of (Pbar_nm cos m lam)^2 over the sphere = 1, so
        // int_{-1}^{1} Pbar_nm Pbar_km dt = (2 - d_m0) * 2 d_nk ... = 2 (m = 0) or 4.
        const int N = 180, Q = 200;
        std::vector<double> x, w;
        gauss_legendre(Q, x, w);
        std::vector<std::vector<double>> tab(Q);
        for (int q = 0; q < Q; ++q) grav::legendre(N, x[q], std::sqrt(1 - x[q] * x[q]), tab[q]);
        double err = 0.0;
        for (int m : {0, 1, 7, 45, 90, 150, 180})
            for (int n = m; n <= N; ++n)
                for (int k = m; k <= N; ++k) {
                    double s = 0.0;
                    for (int q = 0; q < Q; ++q) s += w[q] * tab[q][Field::tri(n, m)] * tab[q][Field::tri(k, m)];
                    double want = (n == k) ? (m == 0 ? 2.0 : 4.0) : 0.0;
                    err = std::max(err, std::fabs(s - want));
                }
        check(err < 1e-11, "orthonormality to degree 180 (200-pt Gauss), max error", err, 0.0);
    }

    std::printf("Acceleration\n");
    {
        std::mt19937_64 rng(68);
        std::uniform_real_distribution<double> U(-1.0, 1.0);
        auto random_point = [&](double rad) {
            Vec3 d{U(rng), U(rng), 0.8 * U(rng)};
            double s = rad / vnorm(d);
            return Vec3{d[0] * s, d[1] * s, d[2] * s};
        };

        double e_sc = 0.0, e_fd = 0.0;
        for (int i = 0; i < 50; ++i) {
            Vec3 x = random_point(egm.R + 200e3 + 40000e3 * (i / 50.0));
            Vec3 a1 = grav::accel_spherical(egm, x), a2 = grav::accel_cunningham(egm, x);
            e_sc = std::max(e_sc, vdiff(a1, a2) / vnorm(a1));
            e_fd = std::max(e_fd, vdiff(a1, fd_gradient(egm, x, 1e-4 * vnorm(x))) / vnorm(a1));
        }
        check(e_sc < 1e-14, "EGM96 4x4: spherical vs Cunningham, max rel", e_sc, 0.0);
        check(e_fd < 1e-10, "EGM96 4x4: spherical vs grad U (4th-order FD), max rel", e_fd, 0.0);

        // Random field to degree 60 with a Kaula-like 1e-5/n^2 spectrum.
        Field rnd(60, egm.GM, egm.R);
        std::normal_distribution<double> G(0.0, 1.0);
        for (int n = 2; n <= 60; ++n)
            for (int m = 0; m <= n; ++m) {
                rnd.c(n, m) = 1e-5 / (n * n) * G(rng);
                rnd.s(n, m) = m ? 1e-5 / (n * n) * G(rng) : 0.0;
            }
        // Drop C00 so both paths and the numerical gradient see only the
        // perturbing part; with the central term in, every comparison would
        // be limited by round-off on a field 1e5 times larger.
        rnd.c(0, 0) = 0.0;
        double e2 = 0.0, e3 = 0.0;
        for (int i = 0; i < 30; ++i) {
            Vec3 x = random_point(rnd.R + 300e3 + 3000e3 * (i / 30.0));
            Vec3 a1 = grav::accel_spherical(rnd, x), a2 = grav::accel_cunningham(rnd, x);
            e2 = std::max(e2, vdiff(a1, a2) / vnorm(a1));
            e3 = std::max(e3, vdiff(a1, fd_gradient(rnd, x, 2e-5 * vnorm(x))) / vnorm(a1));
        }
        check(e2 < 1e-9, "degree-60 perturbing field: spherical vs Cunningham", e2, 0.0);
        check(e3 < 1e-7, "degree-60 perturbing field: spherical vs grad U", e3, 0.0);

        // J2 closed form
        Field j2 = egm.truncated(2);
        j2.c(2, 1) = j2.s(2, 1) = j2.c(2, 2) = j2.s(2, 2) = 0.0;
        double J2 = -std::sqrt(5.0) * j2.c(2, 0), ej = 0.0;
        for (int i = 0; i < 20; ++i) {
            Vec3 x = random_point(j2.R + 500e3 + 1000e3 * i);
            double rr = vnorm(x), k = 1.5 * J2 * std::pow(j2.R / rr, 2), z2 = std::pow(x[2] / rr, 2);
            double mu3 = j2.GM / std::pow(rr, 3);
            Vec3 ref{-mu3 * x[0] * (1 + k * (1 - 5 * z2)), -mu3 * x[1] * (1 + k * (1 - 5 * z2)),
                     -mu3 * x[2] * (1 + k * (3 - 5 * z2))};
            ej = std::max(ej, vdiff(grav::accel_spherical(j2, x), ref) / vnorm(ref));
        }
        check(ej < 1e-14, "J2-only field vs closed form, max rel", ej, 0.0);
    }

    std::printf("GEO equilibria\n");
    {
        // Degree 2 only: a_T = sqrt(15) GM R^2 J22 / r^4 * sin(2 (lam - lam22)) with
        // J22 = sqrt(C22^2 + S22^2); stable on the minor axis lam22 +- 90 deg.
        Field f2 = egm.truncated(2);
        double lam22 = 0.5 * std::atan2(f2.s(2, 2), f2.c(2, 2)) / kDeg;
        double J22 = std::hypot(f2.c(2, 2), f2.s(2, 2));
        auto eq2 = equilibria(f2, r);
        check(eq2.size() == 4, "degree 2: four equilibria", double(eq2.size()), 4.0);
        double pos_err = 0.0, k_err = 0.0;
        double k_ref = 6.0 * std::sqrt(15.0) * f2.GM * f2.R * f2.R * J22 / std::pow(r, 5);
        for (auto& e : eq2) {
            double off = std::remainder(e.lam_deg - lam22, 180.0);
            bool should_be_stable = std::fabs(std::fabs(off) - 90.0) < 1.0;
            double nearest = std::remainder(e.lam_deg - lam22, 90.0);
            pos_err = std::max(pos_err, std::fabs(nearest));
            k_err = std::max(k_err, std::fabs(std::fabs(e.curvature) - k_ref) / k_ref);
            if (e.stable != should_be_stable) ++failures;
        }
        check(pos_err < 1e-6, "degree 2: equilibria at lam22 + k 90 deg (deg)", pos_err, 0.0);
        check(k_err < 1e-6, "degree 2: |d lam_ddot/d lam| = 6 sqrt15 GM R^2 J22 / r^5", k_err, 0.0);

        // Commonly quoted values: stable 75.1 E, 105.3 W; unstable 161.9 E, 11.5 W.
        auto eq = equilibria(egm, r);
        double want[4] = {75.1, 161.9, 360 - 105.3, 360 - 11.5};
        bool stab[4] = {true, false, true, false};
        check(eq.size() == 4, "EGM96 4x4: four equilibria", double(eq.size()), 4.0);
        for (int i = 0; i < 4 && i < int(eq.size()); ++i) {
            char label[80];
            std::snprintf(label, sizeof label, "EGM96 4x4: %s point near %.1f E (within 0.25 deg)",
                          stab[i] ? "stable" : "unstable", want[i]);
            check(std::fabs(eq[i].lam_deg - want[i]) < 0.25 && eq[i].stable == stab[i], label, eq[i].lam_deg, want[i]);
        }
    }

    std::printf("Libration: 1D drift model vs 3D orbit\n");
    {
        auto eq = equilibria(egm, r);
        for (auto& e : eq) {
            if (!e.stable) continue;
            double T0 = 2 * kPi / std::sqrt(-e.curvature);
            double T1 = libration_period_1d(egm, r, (e.lam_deg + 0.5) * kDeg, e.lam_deg * kDeg);
            check_rel("small-angle 2 pi / sqrt(-k) vs 1D at 0.5 deg (days)", T1 / kDay, T0 / kDay, 1e-3);
            double T20 = libration_period_1d(egm, r, (e.lam_deg + 20) * kDeg, e.lam_deg * kDeg);
            Libration3D L = libration_period_3d(egm, (e.lam_deg + 20) * kDeg, e.lam_deg * kDeg, 1);
            check_rel("20 deg release: 3D orbit vs 1D model period (days)", L.period / kDay, T20 / kDay, 1e-3);
        }
    }

    std::printf("East-west station keeping\n");
    {
        double slot = 30.0 * kDeg, half_box = 0.05 * kDeg;
        double pred = std::fabs(along_track(egm, slot, r)) * kYear;
        StationKeeping sk = east_west_keeping(egm, slot, half_box, 2.0);
        check_rel("30 E: simulated dv rate vs |a_T| x 1 yr (m/s/yr)", sk.dv_rate, pred, 0.01);
        double T_pred = 2 * std::sqrt(2 * 2 * half_box / std::fabs(lam_ddot(egm, slot, r))) / kDay;
        check_rel("30 E: cycle period vs 2 sqrt(2 W / A) (days)", sk.cycle_days, T_pred, 0.05);
        check(sk.max_excursion < 0.05, "30 E: daily-mean longitude stays inside +-0.05 deg", sk.max_excursion, 0.05);
    }

    std::printf("\n%s (%d failure%s)\n", failures ? "FAILED" : "ALL PASSED", failures, failures == 1 ? "" : "s");
    return failures ? 1 : 0;
}
