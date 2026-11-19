// Astron Built this AP
// Project 047 | Navigation | GNSS snapshot positioning with RAIM fault detection and exclusion
// Author: Aditya (aditya-astron)
// Date: 2026-11-19 IST
//
// Verification. Run with no arguments for the checks, or with --dump to print
// chi-square quantiles and noncentral CDFs for verify_stats.py to compare
// against scipy.

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <random>

#include "raim.hpp"

using namespace raim;

static int g_fail = 0, g_pass = 0;

static void check(bool ok, const char* what, double got, double want) {
    std::printf("[%s] %-62s got %.6g  want %.6g\n", ok ? "PASS" : "FAIL", what, got, want);
    ok ? ++g_pass : ++g_fail;
}
static void close_rel(const char* what, double got, double want, double rtol) {
    check(std::fabs(got - want) <= rtol * std::fabs(want), what, got, want);
}
static void close_abs(const char* what, double got, double want, double atol) {
    check(std::fabs(got - want) <= atol, what, got, want);
}

static void dump() {
    for (int k : {1, 2, 3, 4, 5, 6, 8, 12})
        for (double p : {0.5, 0.9, 0.99, 0.999, 1.0 - 1e-5, 1.0 - 1e-7})
            std::printf("inv %d %.17g %.17g\n", k, p, chi2_inv(k, p));
    for (int k : {1, 2, 3, 5, 8})
        for (double lam : {0.5, 5.0, 30.0, 120.0})
            for (double x : {1.0, 10.0, 40.0, 150.0}) std::printf("ncx %d %.17g %.17g %.17g\n", k, lam, x, ncchi2_cdf(k, lam, x));
}

// Fixed test sky: the 7-satellite Mumbai epoch from run_raim.
struct Scene {
    std::vector<Vec3> sats;
    std::vector<double> rho;
    Vec3 rx;
    double lat, lon, b;
    Mat H;  // ENU
};

static Scene scene() {
    Scene s;
    s.lat = rad(19.0760);
    s.lon = rad(72.8777);
    s.rx = geodetic_to_ecef(s.lat, s.lon, 14.0);
    s.b = 85e3;
    for (const auto& slot : gps_like_24()) {
        const Vec3 p = sat_ecef(slot, 0.0);
        double az, el;
        az_el(s.rx, s.lat, s.lon, p, az, el);
        if (el < rad(5.0)) continue;
        s.sats.push_back(p);
        s.rho.push_back(true_pseudorange(p, s.rx, s.b));
    }
    s.H = to_enu(solve(s.sats, s.rho).H, s.lat, s.lon);
    return s;
}

static Vec3 enu_err(const Scene& s, const Vec3& r) {
    const auto B = enu_basis(s.lat, s.lon);
    const Vec3 e = sub(r, s.rx);
    Vec3 out;
    for (int k = 0; k < 3; ++k) out[k] = e[0] * B[k][0] + e[1] * B[k][1] + e[2] * B[k][2];
    return out;
}

int main(int argc, char** argv) {
    if (argc > 1 && std::strcmp(argv[1], "--dump") == 0) {
        dump();
        return 0;
    }

    std::printf("-- statistics\n");
    close_rel("chi2_inv(1, 0.95) table value", chi2_inv(1, 0.95), 3.841459, 2e-7);
    close_rel("chi2_inv(5, 0.99) table value", chi2_inv(5, 0.99), 15.086273, 2e-7);
    close_rel("chi2_inv(10, 0.999) table value", chi2_inv(10, 0.999), 29.588298, 2e-7);
    // k = 2: closed form -2 ln(1-p).
    close_rel("chi2_inv(2, 1-1e-7) = -2 ln(1e-7)", chi2_inv(2, 1.0 - 1e-7), -2.0 * std::log(1e-7), 1e-9);
    // Noncentral, 1 dof: P(|Z + sqrt(lam)| < sqrt(x)).
    {
        const double lam = 9.0, x = 16.0;
        const double want = 0.5 * (std::erf((4.0 - 3.0) / std::sqrt(2.0)) + std::erf((4.0 + 3.0) / std::sqrt(2.0)));
        close_rel("ncchi2_cdf(1, 9, 16) vs normal closed form", ncchi2_cdf(1, lam, x), want, 1e-12);
    }
    close_rel("ncchi2_cdf with lambda = 0 equals central", ncchi2_cdf(4, 0.0, 7.0), chi2_cdf(4, 7.0), 1e-15);
    {
        const double t = chi2_inv(3, 1.0 - 1e-5), lam = ncchi2_lambda_for_pmd(3, t, 1e-3);
        close_rel("lambda_det solves P(ncchi2 < T) = Pmd", ncchi2_cdf(3, lam, t), 1e-3, 1e-8);
    }

    std::printf("-- linear algebra and frames\n");
    {
        const Mat A = {{4, 1, 2, 0.5}, {1, 3, 0, 1}, {2, 0, 5, 1}, {0.5, 1, 1, 2}};
        const Mat I = matmul(A, inverse(A));
        double err = 0.0;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) err = std::max(err, std::fabs(I[i][j] - (i == j)));
        close_abs("A * inv(A) = I", err, 0.0, 1e-14);
        // Pole: ECEF z = b = a (1 - f).
        close_rel("geodetic_to_ecef at the pole gives the semi-minor axis",
                  geodetic_to_ecef(kPi / 2, 0, 0)[2], kWgsA * (1.0 - kWgsF), 1e-15);
    }

    std::printf("-- constellation\n");
    {
        const auto c = gps_like_24();
        const double a = 26559.7e3, T = 2.0 * kPi * std::sqrt(a * a * a / kMuE);
        // Nominal a = 26559.7 km sits ~2 km below the exact half-sidereal radius.
        close_abs("period within 10 s of half a sidereal day (s)", T, 0.5 * 86164.0905, 10.0);
        // Same inertial point after one period means ECEF rotated by we T.
        const Vec3 p0 = sat_ecef(c[5], 0.0), p1 = rot_z(sat_ecef(c[5], T), -kOmegaE * T);
        close_abs("circular orbit closes after one period (m)", norm(sub(p0, p1)), 0.0, 1e-4);
        double rmax = 0.0, zmax = 0.0;
        for (const auto& s : c)
            for (double t = 0; t < 43200; t += 600) {
                const Vec3 p = sat_ecef(s, t);
                rmax = std::max(rmax, std::fabs(norm(p) - a));
                zmax = std::max(zmax, p[2] / a);
            }
        close_abs("orbit radius constant (m)", rmax, 0.0, 1e-6);
        close_rel("max |z|/a = sin(55 deg)", zmax, std::sin(rad(55.0)), 2e-4);
    }

    std::printf("-- positioning\n");
    {
        std::mt19937_64 rng(1);
        std::uniform_real_distribution<double> U(0.0, 1.0);
        const auto c = gps_like_24();
        double worst = 0.0;
        int maxit = 0, cases = 0;
        for (int k = 0; k < 300; ++k) {
            const double lat = std::asin(2.0 * U(rng) - 1.0), lon = 2.0 * kPi * U(rng), h = 3000.0 * U(rng);
            const Vec3 r = geodetic_to_ecef(lat, lon, h);
            const double t = 86400.0 * U(rng), b = 3e5 * (U(rng) - 0.5);
            std::vector<Vec3> s;
            std::vector<double> rho;
            for (const auto& sl : c) {
                const Vec3 p = sat_ecef(sl, t);
                double az, el;
                az_el(r, lat, lon, p, az, el);
                if (el < rad(5.0)) continue;
                s.push_back(p);
                rho.push_back(true_pseudorange(p, r, b));
            }
            if (s.size() < 4) continue;
            const Fix f = solve(s, rho);
            worst = std::max(worst, std::max(norm(sub(f.r, r)), std::fabs(f.b - b)));
            maxit = std::max(maxit, f.iterations);
            ++cases;
        }
        close_abs("noiseless fix from Earth's centre, 300 random sites (m)", worst, 0.0, 1e-6);
        check(maxit <= 10, "Gauss-Newton converges within 10 iterations", maxit, 10);

        const Scene sc = scene();
        const Fix nos = solve(sc.sats, sc.rho, false);
        const double sag = norm(sub(nos.r, sc.rx));
        check(sag > 5.0 && sag < 100.0, "skipping Sagnac costs tens of metres", sag, 30.0);

        // ENU rows from the solver equal the az/el construction.
        std::vector<double> az, el;
        for (const auto& p : sc.sats) {
            double a1, e1;
            az_el(sc.rx, sc.lat, sc.lon, p, a1, e1);
            az.push_back(a1);
            el.push_back(e1);
        }
        std::vector<double> rho_ns;
        for (const auto& p : sc.sats) rho_ns.push_back(true_pseudorange(p, sc.rx, sc.b, false));
        const Mat Hae = geometry_enu(az, el), Hs = to_enu(solve(sc.sats, rho_ns, false).H, sc.lat, sc.lon);
        double e = 0.0;
        for (size_t i = 0; i < Hae.size(); ++i)
            for (int k = 0; k < 4; ++k) e = std::max(e, std::fabs(Hae[i][k] - Hs[i][k]));
        close_abs("solver H in ENU matches az/el geometry", e, 0.0, 1e-9);
    }

    std::printf("-- DOP closed form: zenith + ring of n at elevation e\n");
    {
        // H^T H decouples: EE = NN = n cos^2 e / 2, and the (U, clock) block has
        // determinant n (1 - sin e)^2, giving
        // HDOP^2 = 4 / (n cos^2 e), VDOP^2 = (n+1) / (n (1-s)^2), TDOP^2 = (n s^2 + 1) / (n (1-s)^2).
        double worst = 0.0;
        for (int n = 3; n <= 8; ++n)
            for (double ed : {10.0, 25.0, 40.0, 60.0}) {
                std::vector<double> az = {0.0}, el = {kPi / 2};
                for (int k = 0; k < n; ++k) {
                    az.push_back(2.0 * kPi * k / n + 0.3);
                    el.push_back(rad(ed));
                }
                const Dop d = dop(geometry_enu(az, el));
                const double s = std::sin(rad(ed)), c = std::cos(rad(ed));
                worst = std::max({worst, std::fabs(d.hdop / std::sqrt(4.0 / (n * c * c)) - 1.0),
                                  std::fabs(d.vdop / std::sqrt((n + 1.0) / (n * (1 - s) * (1 - s))) - 1.0),
                                  std::fabs(d.tdop / std::sqrt((n * s * s + 1.0) / (n * (1 - s) * (1 - s))) - 1.0)});
            }
        close_abs("HDOP, VDOP, TDOP vs closed form, 24 geometries (rel)", worst, 0.0, 1e-12);
    }

    std::printf("-- RAIM geometry\n");
    const Scene sc = scene();
    {
        const RaimConfig cfg{5.0, 1e-5, 1e-3};
        const RaimGeometry g = raim_geometry(sc.H, cfg);
        double trace = 0.0;
        for (double s : g.s_diag) trace += s;
        close_rel("trace(S) = N - 4", trace, g.dof, 1e-12);
        // Bias B on satellite i moves the noiseless solution by B * A[:, i].
        const int i = g.worst_h;
        const double B = 100.0;
        std::vector<double> r = sc.rho;
        r[i] += B;
        const Vec3 e = enu_err(sc, solve(sc.sats, r).r);
        close_rel("bias response east = B * A_E,i", e[0], B * g.a_east[i], 1e-4);
        close_rel("bias response north = B * A_N,i", e[1], B * g.a_north[i], 1e-4);
        close_rel("bias response up = B * A_U,i", e[2], B * g.a_up[i], 1e-4);
        // And the residual sum of squares is B^2 S_ii.
        const Fix f = solve(sc.sats, r);
        double sse = 0.0;
        for (double x : f.residuals) sse += x * x;
        close_rel("noiseless SSE = B^2 S_ii", sse, B * B * g.s_diag[i], 1e-4);
        // Four satellites: no redundancy, no RAIM.
        const Mat H4(sc.H.begin(), sc.H.begin() + 4);
        check(!raim_geometry(H4, cfg).available, "RAIM unavailable with 4 satellites", 0, 0);
    }

    std::printf("-- Monte Carlo, fault free (100k epochs, sigma 5 m, Pfa 0.01)\n");
    {
        const RaimConfig cfg{5.0, 0.01, 0.05};
        const RaimGeometry g = raim_geometry(sc.H, cfg);
        const Mat Q = inverse(matmul(transpose(sc.H), sc.H));
        std::mt19937_64 rng(2026);
        std::normal_distribution<double> nd(0.0, cfg.sigma);
        const int M = 100000;
        double m[4] = {0, 0, 0, 0}, v[4] = {0, 0, 0, 0}, ts = 0, ts2 = 0;
        int alarms = 0;
        for (int k = 0; k < M; ++k) {
            std::vector<double> r = sc.rho;
            for (double& x : r) x += nd(rng);
            const Fix f = solve(sc.sats, r, true, sc.rx, sc.b);
            const Vec3 e = enu_err(sc, f.r);
            const double ex[4] = {e[0], e[1], e[2], f.b - sc.b};
            for (int j = 0; j < 4; ++j) {
                m[j] += ex[j];
                v[j] += ex[j] * ex[j];
            }
            double sse = 0;
            for (double x : f.residuals) sse += x * x;
            sse /= cfg.sigma * cfg.sigma;
            ts += sse;
            ts2 += sse * sse;
            alarms += sse > g.threshold;
        }
        const char* names[4] = {"var east = sigma^2 Q_EE", "var north = sigma^2 Q_NN", "var up = sigma^2 Q_UU",
                                "var clock = sigma^2 Q_TT"};
        for (int j = 0; j < 4; ++j) {
            const double var = v[j] / M - (m[j] / M) * (m[j] / M);
            close_rel(names[j], var, cfg.sigma * cfg.sigma * Q[j][j], 0.02);
        }
        const double mean = ts / M, var = ts2 / M - mean * mean;
        close_rel("E[SSE/sigma^2] = N - 4", mean, g.dof, 0.01);
        close_rel("Var[SSE/sigma^2] = 2 (N - 4)", var, 2.0 * g.dof, 0.03);
        const double pfa = static_cast<double>(alarms) / M, sd = std::sqrt(0.01 * 0.99 / M);
        close_abs("false-alarm rate = Pfa (5 sigma binomial)", pfa, 0.01, 5 * sd);
    }

    std::printf("-- Monte Carlo, biased satellite (40k epochs each, Pfa 0.01, Pmd 0.05)\n");
    {
        const RaimConfig cfg{5.0, 0.01, 0.05};
        const RaimGeometry g = raim_geometry(sc.H, cfg);
        std::mt19937_64 rng(7);
        std::normal_distribution<double> nd(0.0, cfg.sigma);
        const int M = 40000;
        const double sd = std::sqrt(0.05 * 0.95 / M);
        // At the noncentrality lambda_det every satellite is missed with Pmd,
        // whatever its slope; the worst one also centres its error on the HPL.
        for (int i = 0; i < g.n; ++i) {
            const double B = cfg.sigma * std::sqrt(g.lambda_det / g.s_diag[i]);
            int missed = 0;
            double me = 0, mn = 0;
            for (int k = 0; k < M; ++k) {
                std::vector<double> r = sc.rho;
                for (double& x : r) x += nd(rng);
                r[i] += B;
                const Fix f = solve(sc.sats, r, true, sc.rx, sc.b);
                double sse = 0;
                for (double x : f.residuals) sse += x * x;
                missed += sse / (cfg.sigma * cfg.sigma) <= g.threshold;
                const Vec3 e = enu_err(sc, f.r);
                me += e[0];
                mn += e[1];
            }
            char buf[96];
            std::snprintf(buf, sizeof buf, "missed detection = Pmd, sat %d at lambda_det", i);
            close_abs(buf, static_cast<double>(missed) / M, 0.05, 5 * sd);
            if (i == g.worst_h) close_rel("mean horizontal error at pbias = HPL (worst sat)", std::hypot(me, mn) / M, g.hpl, 0.01);
        }
    }

    std::printf("-- Monte Carlo, integrity: P(herr > HPL and no alarm) <= Pmd over a bias sweep\n");
    {
        const RaimConfig cfg{5.0, 0.01, 0.05};
        const RaimGeometry g = raim_geometry(sc.H, cfg);
        std::mt19937_64 rng(11);
        std::normal_distribution<double> nd(0.0, cfg.sigma);
        const int M = 20000;
        double worst = 0.0;
        for (int i = 0; i < g.n; ++i)
            for (double B = 10.0; B <= 160.0; B += 10.0) {
                int bad = 0;
                for (int k = 0; k < M; ++k) {
                    std::vector<double> r = sc.rho;
                    for (double& x : r) x += nd(rng);
                    r[i] += B;
                    const Fix f = solve(sc.sats, r, true, sc.rx, sc.b);
                    double sse = 0;
                    for (double x : f.residuals) sse += x * x;
                    const Vec3 e = enu_err(sc, f.r);
                    bad += std::hypot(e[0], e[1]) > g.hpl && sse / (cfg.sigma * cfg.sigma) <= g.threshold;
                }
                worst = std::max(worst, static_cast<double>(bad) / M);
            }
        check(worst <= 0.05 + 5 * std::sqrt(0.05 * 0.95 / M), "worst hazardous-misleading rate over 112 cases <= Pmd",
              worst, 0.05);
    }

    std::printf("-- Monte Carlo, exclusion (5k epochs per satellite, 150 m bias, Pfa 1e-5)\n");
    {
        const RaimConfig cfg{5.0, 1e-5, 1e-3};
        std::mt19937_64 rng(13);
        std::normal_distribution<double> nd(0.0, cfg.sigma);
        const int M = 5000;
        int right = 0, wrong = 0, total = 0, clean_alarms = 0;
        for (int i = 0; i < static_cast<int>(sc.sats.size()); ++i)
            for (int k = 0; k < M; ++k) {
                std::vector<double> r = sc.rho;
                for (double& x : r) x += nd(rng);
                const FdeResult clean = detect_and_exclude(sc.sats, r, cfg);
                clean_alarms += clean.alarm;
                r[i] += 150.0;
                const FdeResult fr = detect_and_exclude(sc.sats, r, cfg);
                ++total;
                right += fr.excluded == i;
                wrong += fr.excluded >= 0 && fr.excluded != i;
            }
        close_abs("150 m fault excluded correctly (fraction)", static_cast<double>(right) / total, 1.0, 0.01);
        close_abs("wrong exclusions (fraction)", static_cast<double>(wrong) / total, 0.0, 0.005);
        check(clean_alarms <= 3, "no alarms on the same 35k clean epochs (Pfa 1e-5)", clean_alarms, 0);
    }

    std::printf("\n%d passed, %d failed\n", g_pass, g_fail);
    return g_fail ? 1 : 0;
}
