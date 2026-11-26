// Astron Built this AP
// Project 065 | Algorithmic problem solving / Attitude determination | Lost-in-space star identification with k-vector and Pyramid
// Author: Aditya (aditya-astron)
// Date: 2026-11-26 IST
//
// Lost-in-space star identification for a 12 deg star tracker.
//
//  - Pair catalog: every catalog pair closer than the FOV diagonal, sorted by
//    cos(separation).
//  - k-vector (Mortari 1997): a lookup table over a straight line laid across
//    the sorted values, so a range query costs two index computations and a
//    short trim instead of two binary searches.
//  - Pyramid (Mortari et al. 2004): match a triangle of observed stars with
//    handedness, then confirm it with a fourth star before trusting it.
//  - Davenport q-method for the attitude, Shuster's QUEST covariance for the
//    expected error.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace sid {

constexpr double kPi = 3.14159265358979323846;

struct Vec3 {
    double x = 0, y = 0, z = 0;
};
inline Vec3 operator+(Vec3 a, Vec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3 operator-(Vec3 a, Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3 operator*(double s, Vec3 a) { return {s * a.x, s * a.y, s * a.z}; }
inline double dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(Vec3 a, Vec3 b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm(Vec3 a) { return std::sqrt(dot(a, a)); }
inline Vec3 unit(Vec3 a) { return (1.0 / norm(a)) * a; }
inline double triple(Vec3 a, Vec3 b, Vec3 c) { return dot(a, cross(b, c)); }

// Angle between unit vectors without the acos(dot) precision loss at small angles.
inline double angle(Vec3 a, Vec3 b) { return std::atan2(norm(cross(a, b)), dot(a, b)); }

using Mat3 = std::array<std::array<double, 3>, 3>;

inline Vec3 mul(const Mat3& m, Vec3 v) {
    return {m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
            m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
            m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z};
}
inline Mat3 mul(const Mat3& a, const Mat3& b) {
    Mat3 c{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k) c[i][j] += a[i][k] * b[k][j];
    return c;
}
inline Mat3 transpose(const Mat3& a) {
    Mat3 t{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) t[i][j] = a[j][i];
    return t;
}

// Shuster convention: q = [qv, q4], b = A(q) r.
inline Mat3 quat_to_dcm(const std::array<double, 4>& q) {
    double x = q[0], y = q[1], z = q[2], w = q[3];
    return {{{w * w + x * x - y * y - z * z, 2 * (x * y + w * z), 2 * (x * z - w * y)},
             {2 * (x * y - w * z), w * w - x * x + y * y - z * z, 2 * (y * z + w * x)},
             {2 * (x * z + w * y), 2 * (y * z - w * x), w * w - x * x - y * y + z * z}}};
}

// Small-angle error vector (body frame) of A_est relative to A_true:
// A_est A_true^T = I - [dtheta x] to first order.
inline Vec3 attitude_error(const Mat3& a_est, const Mat3& a_true) {
    Mat3 m = mul(a_est, transpose(a_true));
    return {0.5 * (m[1][2] - m[2][1]), 0.5 * (m[2][0] - m[0][2]), 0.5 * (m[0][1] - m[1][0])};
}

inline Mat3 random_dcm(std::mt19937_64& rng) {
    std::normal_distribution<double> g(0.0, 1.0);
    std::array<double, 4> q{g(rng), g(rng), g(rng), g(rng)};
    double n = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
    for (double& v : q) v /= n;
    return quat_to_dcm(q);
}

// ---------------------------------------------------------------- eigen / q-method

// Cyclic Jacobi for a small symmetric matrix. Returns eigenvalues; columns of v are vectors.
template <int N>
void jacobi_eigen(std::array<std::array<double, N>, N> a, std::array<double, N>& w,
                  std::array<std::array<double, N>, N>& v) {
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) v[i][j] = (i == j);
    for (int sweep = 0; sweep < 60; ++sweep) {
        double off = 0;
        for (int p = 0; p < N; ++p)
            for (int q = p + 1; q < N; ++q) off += a[p][q] * a[p][q];
        if (off < 1e-300) break;
        for (int p = 0; p < N; ++p)
            for (int q = p + 1; q < N; ++q) {
                if (a[p][q] == 0.0) continue;
                double theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
                double t = (theta >= 0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1));
                double c = 1 / std::sqrt(t * t + 1), s = t * c;
                for (int k = 0; k < N; ++k) {
                    double akp = a[k][p], akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for (int k = 0; k < N; ++k) {
                    double apk = a[p][k], aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for (int k = 0; k < N; ++k) {
                    double vkp = v[k][p], vkq = v[k][q];
                    v[k][p] = c * vkp - s * vkq;
                    v[k][q] = s * vkp + c * vkq;
                }
            }
    }
    for (int i = 0; i < N; ++i) w[i] = a[i][i];
}

// Davenport q-method: maximise sum w_i b_i . A r_i. The optimal quaternion is
// the eigenvector of K for its largest eigenvalue.
inline Mat3 q_method(const std::vector<Vec3>& b, const std::vector<Vec3>& r) {
    double B[3][3] = {};
    for (size_t i = 0; i < b.size(); ++i) {
        double bi[3] = {b[i].x, b[i].y, b[i].z}, ri[3] = {r[i].x, r[i].y, r[i].z};
        for (int p = 0; p < 3; ++p)
            for (int q = 0; q < 3; ++q) B[p][q] += bi[p] * ri[q];
    }
    double sig = B[0][0] + B[1][1] + B[2][2];
    double z[3] = {B[1][2] - B[2][1], B[2][0] - B[0][2], B[0][1] - B[1][0]};
    std::array<std::array<double, 4>, 4> K{};
    for (int p = 0; p < 3; ++p) {
        for (int q = 0; q < 3; ++q) K[p][q] = B[p][q] + B[q][p] - (p == q ? sig : 0.0);
        K[p][3] = K[3][p] = z[p];
    }
    K[3][3] = sig;
    std::array<double, 4> w;
    std::array<std::array<double, 4>, 4> v;
    jacobi_eigen<4>(K, w, v);
    int imax = int(std::max_element(w.begin(), w.end()) - w.begin());
    return quat_to_dcm({v[0][imax], v[1][imax], v[2][imax], v[3][imax]});
}

// QUEST covariance for equal-sigma unit vectors (body frame):
// P = sigma^2 [ sum (I - b b^T) ]^-1.  Returns the diagonal.
inline Vec3 quest_cov_diag(const std::vector<Vec3>& b, double sigma) {
    double F[3][3] = {};
    for (const Vec3& v : b) {
        double vv[3] = {v.x, v.y, v.z};
        for (int p = 0; p < 3; ++p)
            for (int q = 0; q < 3; ++q) F[p][q] += (p == q) - vv[p] * vv[q];
    }
    double det = F[0][0] * (F[1][1] * F[2][2] - F[1][2] * F[2][1]) -
                 F[0][1] * (F[1][0] * F[2][2] - F[1][2] * F[2][0]) +
                 F[0][2] * (F[1][0] * F[2][1] - F[1][1] * F[2][0]);
    double s2 = sigma * sigma;
    return {s2 * (F[1][1] * F[2][2] - F[1][2] * F[2][1]) / det,
            s2 * (F[0][0] * F[2][2] - F[0][2] * F[2][0]) / det,
            s2 * (F[0][0] * F[1][1] - F[0][1] * F[1][0]) / det};
}

// ---------------------------------------------------------------- catalog

// Synthetic all-sky catalog, isotropic. 5000 stars is roughly the count to
// visual magnitude 6, which is the regime of a small commercial tracker.
inline std::vector<Vec3> make_catalog(int n, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> g(0.0, 1.0);
    std::vector<Vec3> cat(n);
    for (auto& s : cat) s = unit({g(rng), g(rng), g(rng)});
    return cat;
}

struct Pair {
    double c;  // cos(separation)
    int a, b;
};

// k-vector over a sorted array y[0..n-1]. Line z(j) = q + m j, j = 0..n,
// with z(0) < y_min and z(n) > y_max. K[j] = #{ y_i <= z(j) }.
class KVector {
public:
    KVector() = default;
    explicit KVector(const std::vector<double>& y) : y_(y) {
        n_ = int(y_.size());
        double xi = 4 * std::numeric_limits<double>::epsilon() *
                    std::max(std::fabs(y_.front()), std::fabs(y_.back()));
        q_ = y_.front() - xi;
        m_ = (y_.back() - y_.front() + 2 * xi) / n_;
        K_.assign(n_ + 1, 0);
        int i = 0;
        for (int j = 0; j <= n_; ++j) {
            double zj = q_ + m_ * j;
            while (i < n_ && y_[i] <= zj) ++i;
            K_[j] = i;
        }
    }
    // Half-open index range [first, last) of all y in [ya, yb].
    // `trimmed` reports how many superset elements had to be discarded.
    std::pair<int, int> range(double ya, double yb, int* trimmed = nullptr) const {
        long jlo = long(std::floor((ya - q_) / m_)) - 1;
        long jhi = long(std::ceil((yb - q_) / m_));
        jlo = std::clamp(jlo, 0L, long(n_));
        jhi = std::clamp(jhi, 0L, long(n_));
        int first = K_[jlo], last = K_[jhi];
        int t = 0;
        while (first < last && y_[first] < ya) ++first, ++t;
        while (last > first && y_[last - 1] > yb) --last, ++t;
        if (trimmed) *trimmed = t;
        return {first, last};
    }

private:
    std::vector<double> y_;
    std::vector<int> K_;
    double q_ = 0, m_ = 1;
    int n_ = 0;
};

struct PairCatalog {
    std::vector<Pair> pairs;
    std::vector<double> cosv;
    KVector kv;

    PairCatalog(const std::vector<Vec3>& cat, double max_sep) {
        double cmin = std::cos(max_sep);
        for (int i = 0; i < int(cat.size()); ++i)
            for (int j = i + 1; j < int(cat.size()); ++j) {
                double c = dot(cat[i], cat[j]);
                if (c >= cmin) pairs.push_back({c, i, j});
            }
        std::sort(pairs.begin(), pairs.end(), [](const Pair& p, const Pair& q) { return p.c < q.c; });
        cosv.reserve(pairs.size());
        for (const Pair& p : pairs) cosv.push_back(p.c);
        kv = KVector(cosv);
    }
    // Pairs whose separation lies in [theta - eps, theta + eps].
    std::pair<int, int> query(double theta, double eps, int* trimmed = nullptr) const {
        double lo = std::max(theta - eps, 0.0);
        return kv.range(std::cos(theta + eps), std::cos(lo), trimmed);
    }
};

// ---------------------------------------------------------------- camera

struct Camera {
    int npix = 1024;
    double half_fov = 6.0 * kPi / 180.0;  // square 12 x 12 deg
    double f() const { return 0.5 * npix / std::tan(half_fov); }  // focal length, pixels
    double cx() const { return 0.5 * npix; }
    Vec3 pixel_to_body(double u, double v) const { return unit({(u - cx()) / f(), (v - cx()) / f(), 1.0}); }
    double ifov() const { return 1.0 / f(); }  // rad per pixel at boresight
    // Solid angle of a rectangular pyramid with half-angles a, b: 4 asin(sin a sin b).
    double solid_angle() const { return 4 * std::asin(std::sin(half_fov) * std::sin(half_fov)); }
    double diag_fov() const { return 2 * std::atan(std::sqrt(2.0) * std::tan(half_fov)); }
};

struct Image {
    std::vector<Vec3> b;   // measured body unit vectors
    std::vector<int> truth;  // catalog id or -1 for a spurious detection
};

inline Image render(const std::vector<Vec3>& cat, const Mat3& A, const Camera& cam, double sigma_px,
                    int n_spikes, double drop_prob, std::mt19937_64& rng) {
    std::normal_distribution<double> g(0.0, 1.0);
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    Image im;
    double f = cam.f(), c = cam.cx();
    Vec3 bore{A[2][0], A[2][1], A[2][2]};  // boresight in inertial frame = 3rd row of A
    double cmin = std::cos(cam.diag_fov() / 2 + 1e-3);
    for (int id = 0; id < int(cat.size()); ++id) {
        if (dot(cat[id], bore) < cmin) continue;
        Vec3 v = mul(A, cat[id]);
        double u = f * v.x / v.z + c, w = f * v.y / v.z + c;
        if (u < 0 || u >= cam.npix || w < 0 || w >= cam.npix) continue;
        if (uni(rng) < drop_prob) continue;
        im.b.push_back(cam.pixel_to_body(u + sigma_px * g(rng), w + sigma_px * g(rng)));
        im.truth.push_back(id);
    }
    for (int k = 0; k < n_spikes; ++k) {
        im.b.push_back(cam.pixel_to_body(uni(rng) * cam.npix, uni(rng) * cam.npix));
        im.truth.push_back(-1);
    }
    // Detection order carries no information about identity.
    std::vector<int> idx(im.b.size());
    for (size_t i = 0; i < idx.size(); ++i) idx[i] = int(i);
    std::shuffle(idx.begin(), idx.end(), rng);
    Image out;
    for (int i : idx) out.b.push_back(im.b[i]), out.truth.push_back(im.truth[i]);
    return out;
}

// ---------------------------------------------------------------- Pyramid

struct Solution {
    bool ok = false;
    Mat3 A{};
    std::vector<int> id;  // catalog id per observed star, -1 if unmatched
    int matched = 0;
    long candidates_scanned = 0;
};

class Pyramid {
public:
    Pyramid(const std::vector<Vec3>& cat, const PairCatalog& pc, double eps)
        : cat_(cat), pc_(pc), eps_(eps) {}

    Solution identify(const std::vector<Vec3>& b) const {
        Solution sol;
        int n = int(b.size());
        sol.id.assign(n, -1);
        if (n < 4) return sol;
        // Mortari's index order: walks the spacing (dj, dk) outermost so a
        // single spurious star is abandoned quickly instead of being re-tried
        // in every triangle.
        for (int dj = 1; dj <= n - 2; ++dj)
            for (int dk = 1; dk <= n - dj - 1; ++dk)
                for (int i = 0; i + dj + dk < n; ++i) {
                    int j = i + dj, k = j + dk;
                    std::array<int, 3> tri;
                    int r_used, d_used;
                    if (!match_triangle(b, i, j, k, tri, r_used, d_used, sol.candidates_scanned)) continue;
                    // Coarse attitude from the confirmed four, then label the rest.
                    std::vector<Vec3> bb{b[i], b[j], b[k], b[r_used]};
                    std::vector<Vec3> rr{cat_[tri[0]], cat_[tri[1]], cat_[tri[2]], cat_[d_used]};
                    Mat3 A = q_method(bb, rr);
                    finish(b, A, sol);
                    if (sol.matched >= 4) return sol;
                    sol = Solution{};
                    sol.id.assign(n, -1);
                }
        return sol;
    }

private:
    const std::vector<Vec3>& cat_;
    const PairCatalog& pc_;
    double eps_;

    static uint64_t key(int a, int b) {
        if (a > b) std::swap(a, b);
        return (uint64_t(uint32_t(a)) << 32) | uint32_t(b);
    }

    // Find the catalog triangle (with matching handedness) for observed i,j,k
    // that a fourth observed star confirms. Ambiguity is resolved only by the
    // fourth star; if more than one triangle survives, reject this triple.
    bool match_triangle(const std::vector<Vec3>& b, int i, int j, int k, std::array<int, 3>& out,
                        int& r_used, int& d_used, long& scanned) const {
        double tobs = triple(b[i], b[j], b[k]);
        // Near-collinear triangles cannot be told from their mirror images.
        double sij = angle(b[i], b[j]), sik = angle(b[i], b[k]), sjk = angle(b[j], b[k]);
        if (std::fabs(tobs) < 5 * eps_ * std::max({sij, sik, sjk})) return false;

        auto [a0, a1] = pc_.query(sij, eps_);
        auto [b0, b1] = pc_.query(sik, eps_);
        auto [c0, c1] = pc_.query(sjk, eps_);
        scanned += (a1 - a0) + (b1 - b0) + (c1 - c0);

        std::unordered_map<int, std::vector<int>> nbr_ik;
        for (int p = b0; p < b1; ++p) {
            const Pair& P = pc_.pairs[p];
            nbr_ik[P.a].push_back(P.b);
            nbr_ik[P.b].push_back(P.a);
        }
        std::unordered_set<uint64_t> jk;
        for (int p = c0; p < c1; ++p) jk.insert(key(pc_.pairs[p].a, pc_.pairs[p].b));

        std::vector<std::array<int, 3>> tris;
        for (int p = a0; p < a1; ++p) {
            const Pair& P = pc_.pairs[p];
            for (int o = 0; o < 2; ++o) {
                int ca = o ? P.b : P.a, cb = o ? P.a : P.b;
                auto it = nbr_ik.find(ca);
                if (it == nbr_ik.end()) continue;
                for (int cc : it->second) {
                    if (cc == cb || !jk.count(key(cb, cc))) continue;
                    if ((triple(cat_[ca], cat_[cb], cat_[cc]) > 0) != (tobs > 0)) continue;
                    tris.push_back({ca, cb, cc});
                }
            }
        }
        if (tris.empty()) return false;

        int confirmed = 0;
        for (const auto& t : tris) {
            for (int r = 0; r < int(b.size()); ++r) {
                if (r == i || r == j || r == k) continue;
                int d = confirm(b, i, j, k, r, t);
                if (d >= 0) {
                    if (++confirmed > 1) return false;
                    out = t, r_used = r, d_used = d;
                    break;
                }
            }
        }
        return confirmed == 1;
    }

    // Catalog star d that sits at the observed distances from all three
    // triangle stars, or -1.
    int confirm(const std::vector<Vec3>& b, int i, int j, int k, int r, const std::array<int, 3>& t) const {
        double sir = angle(b[i], b[r]), sjr = angle(b[j], b[r]), skr = angle(b[k], b[r]);
        auto [p0, p1] = pc_.query(sir, eps_);
        int found = -1;
        for (int p = p0; p < p1; ++p) {
            const Pair& P = pc_.pairs[p];
            int d;
            if (P.a == t[0]) d = P.b;
            else if (P.b == t[0]) d = P.a;
            else continue;
            if (d == t[1] || d == t[2]) continue;
            if (std::fabs(angle(cat_[t[1]], cat_[d]) - sjr) > eps_) continue;
            if (std::fabs(angle(cat_[t[2]], cat_[d]) - skr) > eps_) continue;
            if ((triple(cat_[t[0]], cat_[t[1]], cat_[d]) > 0) != (triple(b[i], b[j], b[r]) > 0)) continue;
            if (found >= 0) return -1;  // two fourth stars fit: not a confirmation
            found = d;
        }
        return found;
    }

    // With an attitude in hand every remaining star is a nearest-neighbour
    // lookup; then refine the attitude on all matches.
    void finish(const std::vector<Vec3>& b, const Mat3& A, Solution& sol) const {
        Mat3 At = transpose(A);
        double gate = std::cos(2 * eps_);
        std::vector<Vec3> bb, rr;
        sol.matched = 0;
        for (size_t s = 0; s < b.size(); ++s) {
            Vec3 r = mul(At, b[s]);
            int best = -1;
            double bestc = gate;
            for (int id = 0; id < int(cat_.size()); ++id) {
                double c = dot(r, cat_[id]);
                if (c > bestc) bestc = c, best = id;
            }
            sol.id[s] = best;
        }
        // One catalog star cannot explain two detections. A spike that lands
        // within the gate of a real star claims the same id; drop both rather
        // than guess which one is real.
        std::unordered_map<int, int> uses;
        for (int id : sol.id)
            if (id >= 0) ++uses[id];
        for (size_t s = 0; s < b.size(); ++s) {
            if (sol.id[s] < 0) continue;
            if (uses[sol.id[s]] > 1) {
                sol.id[s] = -1;
                continue;
            }
            bb.push_back(b[s]), rr.push_back(cat_[sol.id[s]]);
            ++sol.matched;
        }
        if (sol.matched >= 4) {
            sol.A = q_method(bb, rr);
            sol.ok = true;
        }
    }
};

}  // namespace sid
