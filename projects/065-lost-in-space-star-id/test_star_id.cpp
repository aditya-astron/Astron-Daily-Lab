// Astron Built this AP
// Project 065 | Algorithmic problem solving / Attitude determination | Lost-in-space star identification with k-vector and Pyramid
// Author: Aditya (aditya-astron)
// Date: 2026-11-26 IST
#include <cstdio>
#include <cstdlib>

#include "star_id.hpp"

using namespace sid;

static int failures = 0;
#define CHECK(cond, ...)                              \
    do {                                              \
        if (!(cond)) {                                \
            std::printf("FAIL %s:%d  ", __FILE__, __LINE__); \
            std::printf(__VA_ARGS__);                 \
            std::printf("\n");                        \
            ++failures;                               \
        }                                             \
    } while (0)

static const double kArcsec = kPi / 180 / 3600;

static void test_kvector_exact() {
    // Duplicates, a dense cluster and both ends of the range.
    std::vector<double> y{-1.0, -0.5, -0.5, -0.5, 0.0, 1e-9, 2e-9, 0.25, 0.25, 0.9, 1.0};
    KVector kv(y);
    std::mt19937_64 rng(3);
    std::uniform_real_distribution<double> u(-1.2, 1.2);
    for (int n = 0; n < 100000; ++n) {
        double a = u(rng), b = u(rng);
        if (a > b) std::swap(a, b);
        if (n % 10 == 0) a = y[n % y.size()];  // hit stored values exactly
        if (n % 7 == 0) b = y[(n / 7) % y.size()];
        if (a > b) std::swap(a, b);
        auto [f, l] = kv.range(a, b);
        long lo = std::lower_bound(y.begin(), y.end(), a) - y.begin();
        long hi = std::upper_bound(y.begin(), y.end(), b) - y.begin();
        CHECK(f == lo && l == hi, "kvector [%g,%g] -> [%d,%d) want [%ld,%ld)", a, b, f, l, lo, hi);
        if (failures > 5) return;
    }
}

static void test_q_method() {
    // Hand case: 90 deg about z, b = A r. A maps x->-y, y->x in Shuster's
    // passive convention: A = [[0,1,0],[-1,0,0],[0,0,1]].
    std::vector<Vec3> r{{1, 0, 0}, {0, 1, 0}}, b{{0, -1, 0}, {1, 0, 0}};
    Mat3 A = q_method(b, r);
    Mat3 want{{{0, 1, 0}, {-1, 0, 0}, {0, 0, 1}}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) CHECK(std::fabs(A[i][j] - want[i][j]) < 1e-12, "q-method hand case");
    // Random rotations, noiseless vectors: exact recovery.
    std::mt19937_64 rng(4);
    auto cat = make_catalog(10, 9);
    for (int n = 0; n < 200; ++n) {
        Mat3 T = random_dcm(rng);
        std::vector<Vec3> bb;
        for (auto& v : cat) bb.push_back(mul(T, v));
        Vec3 e = attitude_error(q_method(bb, cat), T);
        CHECK(norm(e) < 1e-12, "q-method random %d err %g", n, norm(e));
    }
}

static void test_counts(const std::vector<Vec3>& cat, const PairCatalog& pc, const Camera& cam) {
    double N = cat.size();
    double expv = 0.5 * N * (N - 1) * (1 - std::cos(cam.diag_fov())) / 2;
    CHECK(std::fabs(pc.pairs.size() / expv - 1) < 0.02, "pairs %zu vs %g", pc.pairs.size(), expv);
    // Brute force the pair list for a sample of stars.
    for (int i = 0; i < 50; ++i) {
        long brute = 0, listed = 0;
        for (int j = 0; j < int(N); ++j)
            if (j != i && dot(cat[i], cat[j]) >= std::cos(cam.diag_fov())) ++brute;
        for (auto& p : pc.pairs) listed += (p.a == i || p.b == i);
        CHECK(brute == listed, "star %d: %ld pairs listed, %ld brute", i, listed, brute);
    }
    // Solid angle of the 12 x 12 deg pyramid. Closed form: 4 asin(sin^2 6 deg)
    // = 4 asin(0.0109262) = 0.0437057 sr. Independent check: integrate
    // dOmega = dx dy / (1 + x^2 + y^2)^(3/2) over the focal plane, midpoint rule.
    CHECK(std::fabs(cam.solid_angle() - 0.0437057) < 1e-7, "solid angle %.8f", cam.solid_angle());
    double t = std::tan(cam.half_fov), h = 2 * t / 2000, om = 0;
    for (int i = 0; i < 2000; ++i)
        for (int j = 0; j < 2000; ++j) {
            double x = -t + (i + 0.5) * h, y = -t + (j + 0.5) * h;
            om += h * h / std::pow(1 + x * x + y * y, 1.5);
        }
    CHECK(std::fabs(om / cam.solid_angle() - 1) < 1e-7, "solid angle quadrature %.9f", om);
}

static void test_identify(const std::vector<Vec3>& cat, const PairCatalog& pc, const Camera& cam) {
    // Noiseless: every star labelled correctly and attitude exact to round-off.
    {
        Pyramid pyr(cat, pc, 5 * kArcsec);
        std::mt19937_64 rng(5);
        for (int n = 0; n < 300; ++n) {
            Mat3 A = random_dcm(rng);
            Image im = render(cat, A, cam, 0.0, 0, 0.0, rng);
            Solution s = pyr.identify(im.b);
            CHECK(s.ok, "noiseless frame %d not identified (%zu stars)", n, im.b.size());
            if (!s.ok) continue;
            for (size_t i = 0; i < im.b.size(); ++i) CHECK(s.id[i] == im.truth[i], "label frame %d", n);
            CHECK(norm(attitude_error(s.A, A)) < 1e-9, "noiseless attitude");
        }
    }
    // 0.1 px noise with 5 spikes and 20 % drop-outs: never a wrong solution,
    // and RMS error inside 10 % of the QUEST prediction.
    {
        double sig = 0.1 * cam.ifov(), eps = 4 * std::sqrt(2.0) * sig;
        Pyramid pyr(cat, pc, eps);
        std::mt19937_64 rng(6);
        int ok = 0, wrong = 0, M = 1500;
        double sq[3] = {}, pr[3] = {};
        for (int n = 0; n < M; ++n) {
            Mat3 A = random_dcm(rng);
            Image im = render(cat, A, cam, 0.1, 5, 0.2, rng);
            Solution s = pyr.identify(im.b);
            if (!s.ok) continue;
            int right = 0;
            std::vector<Vec3> used;
            for (size_t i = 0; i < im.b.size(); ++i)
                if (s.id[i] >= 0) {
                    used.push_back(im.b[i]);
                    right += s.id[i] == im.truth[i];
                }
            if (right < 4) {
                ++wrong;
                continue;
            }
            ++ok;
            Vec3 e = attitude_error(s.A, A), p = quest_cov_diag(used, sig);
            sq[0] += e.x * e.x, sq[1] += e.y * e.y, sq[2] += e.z * e.z;
            pr[0] += p.x, pr[1] += p.y, pr[2] += p.z;
        }
        CHECK(wrong == 0, "%d wrong identifications", wrong);
        CHECK(ok >= 0.995 * M, "success %d / %d", ok, M);
        for (int k = 0; k < 3; ++k) {
            double ratio = std::sqrt(sq[k] / pr[k]);
            CHECK(std::fabs(ratio - 1) < 0.10, "axis %d RMS / predicted = %.3f", k, ratio);
        }
    }
    // A frame of pure spikes must not produce a solution.
    {
        Pyramid pyr(cat, pc, 4 * std::sqrt(2.0) * 0.1 * cam.ifov());
        std::mt19937_64 rng(8);
        int false_ok = 0;
        for (int n = 0; n < 300; ++n) {
            Image sp = render(cat, random_dcm(rng), cam, 0.1, 12, 1.0, rng);  // every real star dropped
            false_ok += pyr.identify(sp.b).ok;
        }
        CHECK(false_ok == 0, "%d solutions from noise-only frames", false_ok);
    }
}

int main() {
    Camera cam;
    auto cat = make_catalog(5000, 65);
    PairCatalog pc(cat, cam.diag_fov());
    test_kvector_exact();
    test_q_method();
    test_counts(cat, pc, cam);
    test_identify(cat, pc, cam);
    if (failures) {
        std::printf("%d check(s) failed\n", failures);
        return 1;
    }
    std::printf("all checks passed\n");
    return 0;
}
