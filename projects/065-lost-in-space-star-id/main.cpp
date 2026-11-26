// Astron Built this AP
// Project 065 | Algorithmic problem solving / Attitude determination | Lost-in-space star identification with k-vector and Pyramid
// Author: Aditya (aditya-astron)
// Date: 2026-11-26 IST
//
// Study: can a 12 deg, 1024 px tracker with 0.1 px centroids identify the sky
// from nothing, how fast, how often is it wrong, and what accuracy follows.
#include <chrono>
#include <cstdio>

#include "star_id.hpp"

using namespace sid;
using Clock = std::chrono::steady_clock;

static constexpr double kArcsec = kPi / 180.0 / 3600.0;

struct Tally {
    int ok = 0, wrong = 0, none = 0, mislabel = 0;
    double ms = 0;
    double sq[3] = {0, 0, 0}, pred[3] = {0, 0, 0};
};

static Tally run(const std::vector<Vec3>& cat, const PairCatalog& pc, const Camera& cam, double sigma_px,
                 double eps, int spikes, double drop, int trials, uint64_t seed) {
    Pyramid pyr(cat, pc, eps);
    std::mt19937_64 rng(seed);
    Tally t;
    double sig = sigma_px * cam.ifov();
    for (int n = 0; n < trials; ++n) {
        Mat3 A = random_dcm(rng);
        Image im = render(cat, A, cam, sigma_px, spikes, drop, rng);
        auto t0 = Clock::now();
        Solution s = pyr.identify(im.b);
        t.ms += std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
        if (!s.ok) {
            ++t.none;
            continue;
        }
        // Judge the identification by its labels, not by an attitude
        // threshold: with few or clustered stars the roll error of a correct
        // solution can legitimately reach arcminutes.
        int right = 0, bad = 0;
        std::vector<Vec3> used;
        for (size_t i = 0; i < im.b.size(); ++i) {
            if (s.id[i] < 0) continue;
            used.push_back(im.b[i]);
            if (s.id[i] == im.truth[i]) ++right;
            else ++bad;
        }
        if (right < 4) {
            ++t.wrong;
            continue;
        }
        ++t.ok;
        t.mislabel += bad;
        Vec3 e = attitude_error(s.A, A);
        Vec3 p = quest_cov_diag(used, sig);
        t.sq[0] += e.x * e.x, t.sq[1] += e.y * e.y, t.sq[2] += e.z * e.z;
        t.pred[0] += p.x, t.pred[1] += p.y, t.pred[2] += p.z;
    }
    return t;
}

int main() {
    const int N = 5000;
    Camera cam;
    auto cat = make_catalog(N, 65);
    const double sigma_px = 0.1;
    const double sig = sigma_px * cam.ifov();
    const double eps = 4 * std::sqrt(2.0) * sig;  // 4-sigma on a separation

    std::printf("Lost-in-space star ID: %d-star isotropic catalog, %.0f x %.0f deg FOV, %d px\n", N,
                2 * cam.half_fov * 180 / kPi, 2 * cam.half_fov * 180 / kPi, cam.npix);
    std::printf("focal length %.1f px, IFOV %.2f arcsec/px, centroid sigma %.2f px = %.2f arcsec, eps %.1f arcsec\n\n",
                cam.f(), cam.ifov() / kArcsec, sigma_px, sig / kArcsec, eps / kArcsec);

    // 1. Pair catalog size against the isotropic expectation.
    auto t0 = Clock::now();
    PairCatalog pc(cat, cam.diag_fov());
    double build_s = std::chrono::duration<double>(Clock::now() - t0).count();
    double npairs_exp = 0.5 * N * (N - 1) * (1 - std::cos(cam.diag_fov())) / 2;
    std::printf("[1] pair catalog (sep < %.2f deg diagonal)\n", cam.diag_fov() * 180 / kPi);
    std::printf("    pairs %zu, expected N(N-1)/2 (1-cos dmax)/2 = %.0f  (%.2f %%), built in %.2f s, %.1f MB\n",
                pc.pairs.size(), npairs_exp, 100 * (pc.pairs.size() / npairs_exp - 1), build_s,
                pc.pairs.size() * (sizeof(Pair) + sizeof(double) + sizeof(int)) / 1e6);

    // 2. Stars per frame against N * Omega / 4 pi.
    {
        std::mt19937_64 rng(1);
        double sum = 0;
        int below4 = 0, M = 20000;
        for (int n = 0; n < M; ++n) {
            Image im = render(cat, random_dcm(rng), cam, 0.0, 0, 0.0, rng);
            sum += im.b.size();
            below4 += im.b.size() < 4;
        }
        double lam = N * cam.solid_angle() / (4 * kPi);
        double p4 = std::exp(-lam) * (1 + lam + lam * lam / 2 + lam * lam * lam / 6);
        std::printf("\n[2] stars in frame: mean %.2f (N Omega/4pi = %.2f), P(<4) %.4f (Poisson %.4f)\n",
                    sum / M, lam, double(below4) / M, p4);
    }

    // 3. k-vector range search: correctness, candidate count, speed.
    {
        std::mt19937_64 rng(2);
        std::uniform_real_distribution<double> uth(0.2 * kPi / 180, cam.diag_fov() - eps);
        const int Q = 200000;
        long trimmed = 0, mismatch = 0;
        std::vector<std::pair<double, double>> qs(Q);
        for (auto& q : qs) q = {uth(rng), eps};
        for (auto [th, e] : qs) {
            int tr;
            auto [a, b] = pc.query(th, e, &tr);
            trimmed += tr;
            auto lo = std::lower_bound(pc.cosv.begin(), pc.cosv.end(), std::cos(th + e)) - pc.cosv.begin();
            auto hi = std::upper_bound(pc.cosv.begin(), pc.cosv.end(), std::cos(th - e)) - pc.cosv.begin();
            mismatch += (a != lo) || (b != hi);
        }
        long sink = 0;
        auto t1 = Clock::now();
        for (auto [th, e] : qs) sink += pc.query(th, e).first;
        double kv_ns = std::chrono::duration<double, std::nano>(Clock::now() - t1).count() / Q;
        t1 = Clock::now();
        for (auto [th, e] : qs)
            sink += std::lower_bound(pc.cosv.begin(), pc.cosv.end(), std::cos(th + e)) - pc.cosv.begin() +
                    (std::upper_bound(pc.cosv.begin(), pc.cosv.end(), std::cos(th - e)) - pc.cosv.begin());
        double bs_ns = std::chrono::duration<double, std::nano>(Clock::now() - t1).count() / Q;
        std::printf("\n[3] k-vector vs binary search over %d random queries\n", Q);
        std::printf("    index mismatches %ld, superset elements trimmed per query %.2f\n", mismatch,
                    double(trimmed) / Q);
        std::printf("    %.0f ns per k-vector query vs %.0f ns for two binary searches (sink %ld)\n", kv_ns, bs_ns,
                    sink % 7);
        // Pair separations on an isotropic sphere have pdf sin(t)/2, so a
        // +-eps window at t should hold N(N-1)/2 (cos(t-e) - cos(t+e))/2 pairs.
        std::printf("    candidates per side query, mean over 20000 queries per band vs isotropic expectation:\n");
        for (auto [lo, hi] : std::vector<std::pair<double, double>>{{0.5, 2}, {2, 5}, {5, 8}, {8, 12}, {12, 16.5}}) {
            std::uniform_real_distribution<double> ub(lo * kPi / 180, hi * kPi / 180);
            double got = 0, expv = 0;
            for (int n = 0; n < 20000; ++n) {
                double th = ub(rng);
                auto [a, b] = pc.query(th, eps);
                got += b - a;
                expv += 0.5 * N * (N - 1) * (std::cos(th - eps) - std::cos(th + eps)) / 2;
            }
            std::printf("      %4.1f-%4.1f deg: %6.1f  expected %6.1f  (%+.2f %%)\n", lo, hi, got / 20000,
                        expv / 20000, 100 * (got / expv - 1));
        }
    }

    // 4. Robustness: spurious detections and dropped stars.
    std::printf("\n[4] Pyramid, %.2f px noise, 2000 frames per row\n", sigma_px);
    std::printf("    spikes  drop   correct    wrong  no-ID  mislabelled   ms/frame\n");
    for (auto [sp, dr] : std::vector<std::pair<int, double>>{{0, 0}, {2, 0}, {5, 0}, {10, 0}, {5, 0.2}, {10, 0.3}}) {
        Tally t = run(cat, pc, cam, sigma_px, eps, sp, dr, 2000, 100 + sp);
        std::printf("    %6d  %4.1f  %8.2f%%  %6d  %5d  %11d  %9.2f\n", sp, dr, 100.0 * t.ok / 2000, t.wrong,
                    t.none, t.mislabel, t.ms / 2000);
    }

    // 5. Accuracy against Shuster's QUEST covariance.
    std::printf("\n[5] attitude error vs QUEST covariance (3000 frames, no spikes)\n");
    std::printf("    sigma_px   RMS x / y / z (arcsec)        predicted x / y / z       correct  no-ID\n");
    for (double spx : {0.05, 0.1, 0.3, 1.0}) {
        double e = 4 * std::sqrt(2.0) * spx * cam.ifov();
        Tally t = run(cat, pc, cam, spx, e, 0, 0.0, 3000, 7);
        double m = t.ok;
        std::printf("    %7.2f   %6.2f %6.2f %7.2f        %6.2f %6.2f %7.2f      %6.2f%%  %5d\n", spx,
                    std::sqrt(t.sq[0] / m) / kArcsec, std::sqrt(t.sq[1] / m) / kArcsec,
                    std::sqrt(t.sq[2] / m) / kArcsec, std::sqrt(t.pred[0] / m) / kArcsec,
                    std::sqrt(t.pred[1] / m) / kArcsec, std::sqrt(t.pred[2] / m) / kArcsec, 100.0 * m / 3000, t.none);
    }
    return 0;
}
