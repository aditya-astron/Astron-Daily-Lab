// Astron Built this AP
// Project 047 | Navigation | GNSS snapshot positioning with RAIM fault detection and exclusion
// Author: Aditya (aditya-astron)
// Date: 2026-11-19 IST

#include "raim.hpp"

#include <algorithm>
#include <limits>
#include <map>
#include <tuple>

namespace raim {

// ---------------------------------------------------------------- frames

Vec3 geodetic_to_ecef(double lat, double lon, double h) {
    const double e2 = kWgsF * (2.0 - kWgsF);
    const double sl = std::sin(lat), cl = std::cos(lat);
    const double N = kWgsA / std::sqrt(1.0 - e2 * sl * sl);
    return {(N + h) * cl * std::cos(lon), (N + h) * cl * std::sin(lon), (N * (1.0 - e2) + h) * sl};
}

std::array<Vec3, 3> enu_basis(double lat, double lon) {
    const double sl = std::sin(lat), cl = std::cos(lat), so = std::sin(lon), co = std::cos(lon);
    return {Vec3{-so, co, 0.0}, Vec3{-sl * co, -sl * so, cl}, Vec3{cl * co, cl * so, sl}};
}

Vec3 rot_z(const Vec3& v, double a) {
    const double c = std::cos(a), s = std::sin(a);
    return {c * v[0] + s * v[1], -s * v[0] + c * v[1], v[2]};
}

// ------------------------------------------------------------ linear algebra

Mat matmul(const Mat& a, const Mat& b) {
    const size_t n = a.size(), m = b[0].size(), k = b.size();
    Mat c(n, std::vector<double>(m, 0.0));
    for (size_t i = 0; i < n; ++i)
        for (size_t p = 0; p < k; ++p) {
            const double aip = a[i][p];
            for (size_t j = 0; j < m; ++j) c[i][j] += aip * b[p][j];
        }
    return c;
}

Mat transpose(const Mat& a) {
    Mat t(a[0].size(), std::vector<double>(a.size()));
    for (size_t i = 0; i < a.size(); ++i)
        for (size_t j = 0; j < a[0].size(); ++j) t[j][i] = a[i][j];
    return t;
}

Mat inverse(Mat a) {
    const size_t n = a.size();
    Mat inv(n, std::vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i) inv[i][i] = 1.0;
    for (size_t c = 0; c < n; ++c) {
        size_t p = c;
        for (size_t r = c + 1; r < n; ++r)
            if (std::fabs(a[r][c]) > std::fabs(a[p][c])) p = r;
        if (std::fabs(a[p][c]) < 1e-300) throw std::runtime_error("singular matrix");
        std::swap(a[p], a[c]);
        std::swap(inv[p], inv[c]);
        const double d = a[c][c];
        for (size_t j = 0; j < n; ++j) {
            a[c][j] /= d;
            inv[c][j] /= d;
        }
        for (size_t r = 0; r < n; ++r) {
            if (r == c) continue;
            const double f = a[r][c];
            if (f == 0.0) continue;
            for (size_t j = 0; j < n; ++j) {
                a[r][j] -= f * a[c][j];
                inv[r][j] -= f * inv[c][j];
            }
        }
    }
    return inv;
}

// ------------------------------------------------------------- statistics

// P(a,x) and Q(a,x): series below x = a+1, Lentz continued fraction above, so
// whichever tail is small is computed directly rather than as 1 - (near 1).
static void gamma_pq(double a, double x, double& P, double& Q) {
    if (x <= 0.0) {
        P = 0.0;
        Q = 1.0;
        return;
    }
    const double lnpre = -x + a * std::log(x) - std::lgamma(a);
    if (x < a + 1.0) {
        double ap = a, del = 1.0 / a, sum = del;
        for (int n = 0; n < 10000; ++n) {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if (std::fabs(del) < std::fabs(sum) * 1e-17) break;
        }
        P = sum * std::exp(lnpre);
        Q = 1.0 - P;
    } else {
        const double tiny = 1e-300;
        double b = x + 1.0 - a, c = 1.0 / tiny, d = 1.0 / b, h = d;
        for (int i = 1; i < 10000; ++i) {
            const double an = -i * (i - a);
            b += 2.0;
            d = an * d + b;
            if (std::fabs(d) < tiny) d = tiny;
            c = b + an / c;
            if (std::fabs(c) < tiny) c = tiny;
            d = 1.0 / d;
            const double del = d * c;
            h *= del;
            if (std::fabs(del - 1.0) < 1e-17) break;
        }
        Q = std::exp(lnpre) * h;
        P = 1.0 - Q;
    }
}

double gamma_p(double a, double x) {
    double P, Q;
    gamma_pq(a, x, P, Q);
    return P;
}

double chi2_cdf(double k, double x) { return gamma_p(0.5 * k, 0.5 * x); }

double chi2_inv(double k, double p) {
    if (p <= 0.0) return 0.0;
    if (p >= 1.0) return INFINITY;
    // Bisect on whichever tail is small to keep relative accuracy near p -> 1.
    const bool upper = p > 0.5;
    const double target = upper ? 1.0 - p : p;
    auto f = [&](double x) {
        double P, Q;
        gamma_pq(0.5 * k, 0.5 * x, P, Q);
        return upper ? target - Q : P - target;  // increasing in x
    };
    double lo = 0.0, hi = std::max(1.0, k);
    while (f(hi) < 0.0) hi *= 2.0;
    for (int i = 0; i < 200 && hi - lo > 1e-14 * hi; ++i) {
        const double mid = 0.5 * (lo + hi);
        (f(mid) < 0.0 ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
}

// Poisson mixture of central chi-squares: sum_j Pois(j; lambda/2) P(k/2 + j, x/2),
// summed outward from the Poisson mode until the weights are negligible.
double ncchi2_cdf(double k, double lambda, double x) {
    if (lambda <= 0.0) return chi2_cdf(k, x);
    const double m = 0.5 * lambda;
    const int j0 = static_cast<int>(std::floor(m));
    auto w = [&](int j) { return std::exp(-m + j * std::log(m) - std::lgamma(j + 1.0)); };
    double sum = 0.0;
    for (int j = j0;; ++j) {
        const double wj = w(j);
        sum += wj * gamma_p(0.5 * k + j, 0.5 * x);
        if (j > j0 + 5 && wj < 1e-18) break;
    }
    for (int j = j0 - 1; j >= 0; --j) {
        const double wj = w(j);
        sum += wj * gamma_p(0.5 * k + j, 0.5 * x);
        if (wj < 1e-18) break;
    }
    return std::min(1.0, sum);
}

double ncchi2_lambda_for_pmd(double k, double t, double pmd) {
    if (chi2_cdf(k, t) <= pmd) return 0.0;
    double lo = 0.0, hi = 1.0;
    while (ncchi2_cdf(k, hi, t) > pmd) hi *= 2.0;
    for (int i = 0; i < 200 && hi - lo > 1e-12 * hi; ++i) {
        const double mid = 0.5 * (lo + hi);
        (ncchi2_cdf(k, mid, t) > pmd ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
}

// ------------------------------------------------------------- positioning

// Satellite position rotated into the receive-time ECEF frame. tau is
// iterated because it depends on the rotated position itself.
static Vec3 sagnac_sat(const Vec3& sat, const Vec3& r) {
    double tau = norm(sub(sat, r)) / kC;
    Vec3 s = sat;
    for (int i = 0; i < 3; ++i) {
        s = rot_z(sat, kOmegaE * tau);
        tau = norm(sub(s, r)) / kC;
    }
    return s;
}

double true_pseudorange(const Vec3& sat, const Vec3& r, double b, bool sagnac) {
    const Vec3 s = sagnac ? sagnac_sat(sat, r) : sat;
    return norm(sub(s, r)) + b;
}

Fix solve(const std::vector<Vec3>& sats, const std::vector<double>& rho, bool sagnac, Vec3 r0,
          double b0, double tol, int max_iter) {
    const size_t n = sats.size();
    if (n < 4 || rho.size() != n) throw std::invalid_argument("need >= 4 matched pseudoranges");
    Fix f;
    f.r = r0;
    f.b = b0;
    f.H.assign(n, std::vector<double>(4));
    f.residuals.assign(n, 0.0);
    for (int it = 1; it <= max_iter; ++it) {
        for (size_t i = 0; i < n; ++i) {
            const Vec3 s = sagnac ? sagnac_sat(sats[i], f.r) : sats[i];
            const Vec3 d = sub(s, f.r);
            const double rg = norm(d);
            f.H[i] = {-d[0] / rg, -d[1] / rg, -d[2] / rg, 1.0};
            f.residuals[i] = rho[i] - (rg + f.b);
        }
        // Normal equations; H^T H is 4x4 and well conditioned for any usable sky.
        const Mat Ht = transpose(f.H);
        const Mat N = inverse(matmul(Ht, f.H));
        Mat y(n, std::vector<double>(1));
        for (size_t i = 0; i < n; ++i) y[i][0] = f.residuals[i];
        const Mat dx = matmul(N, matmul(Ht, y));
        for (int k = 0; k < 3; ++k) f.r[k] += dx[k][0];
        f.b += dx[3][0];
        f.iterations = it;
        const double step = std::sqrt(dx[0][0] * dx[0][0] + dx[1][0] * dx[1][0] + dx[2][0] * dx[2][0] +
                                      dx[3][0] * dx[3][0]);
        if (step < tol) break;
    }
    // Final residuals and geometry at the converged point.
    for (size_t i = 0; i < n; ++i) {
        const Vec3 s = sagnac ? sagnac_sat(sats[i], f.r) : sats[i];
        const Vec3 d = sub(s, f.r);
        const double rg = norm(d);
        f.H[i] = {-d[0] / rg, -d[1] / rg, -d[2] / rg, 1.0};
        f.residuals[i] = rho[i] - (rg + f.b);
    }
    return f;
}

Mat geometry_enu(const std::vector<double>& az, const std::vector<double>& el) {
    Mat H(az.size(), std::vector<double>(4));
    for (size_t i = 0; i < az.size(); ++i)
        H[i] = {-std::cos(el[i]) * std::sin(az[i]), -std::cos(el[i]) * std::cos(az[i]), -std::sin(el[i]), 1.0};
    return H;
}

Mat to_enu(const Mat& H_ecef, double lat, double lon) {
    const auto B = enu_basis(lat, lon);
    Mat H(H_ecef.size(), std::vector<double>(4));
    for (size_t i = 0; i < H_ecef.size(); ++i) {
        for (int k = 0; k < 3; ++k)
            H[i][k] = H_ecef[i][0] * B[k][0] + H_ecef[i][1] * B[k][1] + H_ecef[i][2] * B[k][2];
        H[i][3] = H_ecef[i][3];
    }
    return H;
}

Dop dop(const Mat& H_enu) {
    const Mat Q = inverse(matmul(transpose(H_enu), H_enu));
    Dop d;
    d.hdop = std::sqrt(Q[0][0] + Q[1][1]);
    d.vdop = std::sqrt(Q[2][2]);
    d.tdop = std::sqrt(Q[3][3]);
    d.pdop = std::sqrt(Q[0][0] + Q[1][1] + Q[2][2]);
    d.gdop = std::sqrt(Q[0][0] + Q[1][1] + Q[2][2] + Q[3][3]);
    return d;
}

// --------------------------------------------------------------- integrity

RaimGeometry raim_geometry(const Mat& H, const RaimConfig& cfg) {
    RaimGeometry g;
    g.n = static_cast<int>(H.size());
    g.dof = g.n - 4;
    if (g.dof < 1) return g;
    const Mat A = matmul(inverse(matmul(transpose(H), H)), transpose(H));  // 4 x N
    // Threshold and lambda depend only on (dof, Pfa, Pmd); availability sweeps
    // call this thousands of times, so solve each combination once.
    static std::map<std::tuple<int, double, double>, std::pair<double, double>> cache;
    const auto key = std::make_tuple(g.dof, cfg.pfa, cfg.pmd);
    auto it = cache.find(key);
    if (it == cache.end()) {
        const double t = chi2_inv(g.dof, 1.0 - cfg.pfa);
        it = cache.emplace(key, std::make_pair(t, ncchi2_lambda_for_pmd(g.dof, t, cfg.pmd))).first;
    }
    g.threshold = it->second.first;
    g.lambda_det = it->second.second;
    double hmax = 0.0, vmax = 0.0;
    for (int i = 0; i < g.n; ++i) {
        double hA = 0.0;  // (H A)_ii
        for (int k = 0; k < 4; ++k) hA += H[i][k] * A[k][i];
        const double s = 1.0 - hA;
        g.s_diag.push_back(s);
        g.a_east.push_back(A[0][i]);
        g.a_north.push_back(A[1][i]);
        g.a_up.push_back(A[2][i]);
        // A satellite with S_ii -> 0 cannot be tested: its bias never shows in the residuals.
        const double hs = s > 1e-12 ? std::hypot(A[0][i], A[1][i]) / std::sqrt(s) : INFINITY;
        const double vs = s > 1e-12 ? std::fabs(A[2][i]) / std::sqrt(s) : INFINITY;
        g.hslope.push_back(hs);
        g.vslope.push_back(vs);
        if (hs > hmax) {
            hmax = hs;
            g.worst_h = i;
        }
        vmax = std::max(vmax, vs);
    }
    const double pbias = cfg.sigma * std::sqrt(g.lambda_det);
    g.hpl = hmax * pbias;
    g.vpl = vmax * pbias;
    g.available = true;
    return g;
}

static double sse_over_sigma2(const Fix& f, double sigma) {
    double s = 0.0;
    for (double r : f.residuals) s += r * r;
    return s / (sigma * sigma);
}

FdeResult detect_and_exclude(const std::vector<Vec3>& sats, const std::vector<double>& rho,
                             const RaimConfig& cfg, bool sagnac) {
    FdeResult out;
    const size_t n = sats.size();
    if (n < 5) throw std::invalid_argument("detection needs >= 5 satellites");
    const Fix full = solve(sats, rho, sagnac);
    out.test = sse_over_sigma2(full, cfg.sigma);
    out.alarm = out.test > chi2_inv(static_cast<double>(n) - 4.0, 1.0 - cfg.pfa);
    if (!out.alarm || n < 6) return out;
    const double t_sub = chi2_inv(static_cast<double>(n) - 5.0, 1.0 - cfg.pfa);
    double best = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < n; ++i) {
        std::vector<Vec3> s2;
        std::vector<double> r2;
        for (size_t j = 0; j < n; ++j)
            if (j != i) {
                s2.push_back(sats[j]);
                r2.push_back(rho[j]);
            }
        const Fix fi = solve(s2, r2, sagnac, full.r, full.b);
        const double ti = sse_over_sigma2(fi, cfg.sigma);
        if (ti < t_sub && ti < best) {
            best = ti;
            out.excluded = static_cast<int>(i);
        }
    }
    out.recovered = out.excluded >= 0;
    return out;
}

// --------------------------------------------------------------- constellation

std::vector<Slot> gps_like_24() {
    // Plane RAANs 60 deg apart; uneven in-plane slots so that a single
    // missing satellite leaves a real hole, as in the operational layout.
    const double raan[6] = {272.847, 332.847, 32.847, 92.847, 152.847, 212.847};
    const double u[6][4] = {{268.126, 161.786, 11.676, 41.806},   {80.956, 173.336, 309.976, 204.376},
                            {111.876, 11.796, 339.666, 241.556},  {135.226, 265.446, 35.156, 167.356},
                            {197.046, 302.596, 66.066, 333.686},  {238.886, 345.226, 105.206, 135.346}};
    std::vector<Slot> s;
    for (int p = 0; p < 6; ++p)
        for (int k = 0; k < 4; ++k) s.push_back({rad(raan[p]), rad(u[p][k])});
    return s;
}

Vec3 sat_ecef(const Slot& s, double t, double incl, double a) {
    const double n = std::sqrt(kMuE / (a * a * a));
    const double u = s.u0 + n * t;
    const double cO = std::cos(s.raan), sO = std::sin(s.raan), ci = std::cos(incl), si = std::sin(incl);
    const Vec3 eci = {a * (cO * std::cos(u) - sO * std::sin(u) * ci), a * (sO * std::cos(u) + cO * std::sin(u) * ci),
                      a * std::sin(u) * si};
    return rot_z(eci, kOmegaE * t);
}

void az_el(const Vec3& r_rx, double lat, double lon, const Vec3& r_sat, double& az, double& el) {
    const auto B = enu_basis(lat, lon);
    const Vec3 d = sub(r_sat, r_rx);
    const double rg = norm(d);
    const double e = (d[0] * B[0][0] + d[1] * B[0][1] + d[2] * B[0][2]) / rg;
    const double nn = (d[0] * B[1][0] + d[1] * B[1][1] + d[2] * B[1][2]) / rg;
    const double u = (d[0] * B[2][0] + d[1] * B[2][1] + d[2] * B[2][2]) / rg;
    el = std::asin(u);
    az = std::atan2(e, nn);
    if (az < 0.0) az += 2.0 * kPi;
}

}  // namespace raim
