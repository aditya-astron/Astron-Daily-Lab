// Astron Built this AP
// Project 021 | Algorithmic problem solving | Ground-station antenna scheduling by min-cost flow
// Author: Aditya (aditya-astron)
// Date: 2026-11-07 IST
//
// Verification:
//   1. Pass geometry against hand-computed overhead-pass durations.
//   2. A hand-worked instance where greedy is provably suboptimal.
//   3. Flow vs. exhaustive subset search on 600 random instances, k = 1..3.
//   4. Flow vs. weighted-interval DP for k = 1 on larger instances.
//   5. Every returned schedule is feasible: no antenna double-booked within
//      the turnaround, at most k antennas used, value matches the passes.
//   6. Optimal value is monotone in k and saturates at the peak concurrency.

#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "scheduler.hpp"

using sched::Pass;

static int failures = 0;
#define CHECK(cond, ...)                                       \
    do {                                                       \
        if (!(cond)) {                                         \
            ++failures;                                        \
            std::printf("FAIL %s:%d: ", __FILE__, __LINE__);   \
            std::printf(__VA_ARGS__);                          \
            std::printf("\n");                                 \
        }                                                      \
    } while (0)

static int64_t brute_force(const std::vector<Pass>& p, int k, int64_t tau) {
    const int n = static_cast<int>(p.size());
    int64_t best = 0;
    for (uint32_t m = 0; m < (1u << n); ++m) {
        std::vector<int> sub;
        int64_t v = 0;
        for (int j = 0; j < n; ++j)
            if (m >> j & 1) {
                sub.push_back(j);
                v += p[j].value;
            }
        if (v > best && sched::max_overlap(p, sub, tau) <= k) best = v;
    }
    return best;
}

static bool feasible(const std::vector<Pass>& p, const sched::Schedule& s, int k, int64_t tau) {
    int64_t v = 0;
    for (size_t i = 0; i < p.size(); ++i) {
        if (s.antenna[i] < -1 || s.antenna[i] >= k) return false;
        if (s.antenna[i] < 0) continue;
        v += p[i].value;
        for (size_t j = i + 1; j < p.size(); ++j) {
            if (s.antenna[j] != s.antenna[i]) continue;
            bool overlap = p[i].aos < sched::busy_end(p[j], tau) && p[j].aos < sched::busy_end(p[i], tau);
            if (overlap) return false;
        }
    }
    return v == s.value;
}

static std::vector<Pass> random_passes(std::mt19937_64& rng, int n, int64_t horizon) {
    std::uniform_int_distribution<int64_t> start(0, horizon), dur(60, 900), val(1, 1000);
    std::vector<Pass> p(n);
    for (int j = 0; j < n; ++j) {
        p[j].sat = j;
        p[j].aos = start(rng);
        p[j].los = p[j].aos + dur(rng);
        p[j].value = val(rng);
    }
    return p;
}

int main() {
    // 1. Geometry. Overhead pass, 0 deg mask, 500 km: half-arc = acos(Re/r),
    //    duration = 2 acos(6378.137/6878.137) / n = 693.3 s (computed by hand
    //    with n = sqrt(mu/r^3) = 1.10678e-3 rad/s).
    {
        auto g = sched::pass_geometry(500, 0.0, 0.0);
        CHECK(std::fabs(g.duration_s - 693.3) < 0.5, "overhead 0-deg pass %.2f s", g.duration_s);
        CHECK(std::fabs(g.max_elev_rad - M_PI / 2) < 1e-12, "zenith elevation");
        // 10 deg mask shortens it: lam = acos(Re cos10 / r) - 10 deg = 0.24515 rad
        // -> 2 * 0.24515 / 1.10678e-3 = 443.0 s.
        auto g10 = sched::pass_geometry(500, 10 * M_PI / 180, 0.0);
        CHECK(std::fabs(g10.duration_s - 443.0) < 0.5, "overhead 10-deg pass %.2f s", g10.duration_s);
        // A grazing pass at the mask limit has zero length and peaks at the mask.
        double lam = sched::earth_central_angle_limit(500, 10 * M_PI / 180);
        auto gg = sched::pass_geometry(500, 10 * M_PI / 180, lam * (1 - 1e-9));
        CHECK(gg.duration_s < 1.0, "grazing pass %.3f s", gg.duration_s);
        CHECK(std::fabs(gg.max_elev_rad - 10 * M_PI / 180) < 1e-6, "grazing peak elev %.6f", gg.max_elev_rad);
    }

    // 2. Hand case, one antenna, no turnaround. A long valuable pass (10) sits
    //    across two shorter ones (7 + 6). Greedy grabs the 10; optimum is 13.
    {
        std::vector<Pass> p = {{0, 0, 100, 10}, {1, 0, 50, 7}, {2, 50, 100, 6}};
        auto ex = sched::schedule_mincost_flow(p, 1, 0);
        auto gr = sched::schedule_greedy(p, 1, 0);
        CHECK(ex.value == 13, "hand case exact %lld", (long long)ex.value);
        CHECK(gr.value == 10, "hand case greedy %lld", (long long)gr.value);
        // With 60 s turnaround, passes 1 and 2 now collide, so the best is 10.
        auto ex_t = sched::schedule_mincost_flow(p, 1, 60);
        CHECK(ex_t.value == 10, "hand case with turnaround %lld", (long long)ex_t.value);
        // Two antennas: take 10 on one, and the better of 7/6 chain on the other (both fit).
        auto ex2 = sched::schedule_mincost_flow(p, 2, 0);
        CHECK(ex2.value == 23, "hand case k=2 %lld", (long long)ex2.value);
    }

    // 3. Exhaustive search on small random instances.
    std::mt19937_64 rng(7);
    int greedy_suboptimal = 0, cases = 0;
    for (int trial = 0; trial < 200; ++trial) {
        int n = 6 + trial % 9;  // 6..14 passes
        auto p = random_passes(rng, n, 3600);
        int64_t tau = (trial % 3) * 60;
        for (int k = 1; k <= 3; ++k) {
            auto ex = sched::schedule_mincost_flow(p, k, tau);
            auto gr = sched::schedule_greedy(p, k, tau);
            int64_t bf = brute_force(p, k, tau);
            ++cases;
            CHECK(ex.value == bf, "trial %d k=%d flow %lld brute %lld", trial, k, (long long)ex.value,
                  (long long)bf);
            CHECK(feasible(p, ex, k, tau), "trial %d k=%d flow schedule infeasible", trial, k);
            CHECK(feasible(p, gr, k, tau), "trial %d k=%d greedy schedule infeasible", trial, k);
            CHECK(gr.value <= ex.value, "greedy beat the optimum?");
            greedy_suboptimal += gr.value < ex.value;
        }
    }
    std::printf("exhaustive: %d cases match, greedy suboptimal in %d\n", cases, greedy_suboptimal);

    // 4. DP cross-check for k = 1 at sizes brute force cannot reach.
    for (int trial = 0; trial < 50; ++trial) {
        auto p = random_passes(rng, 300, 86400);
        auto ex = sched::schedule_mincost_flow(p, 1, 120);
        auto dp = sched::schedule_dp_single(p, 120);
        CHECK(ex.value == dp.value, "dp trial %d flow %lld dp %lld", trial, (long long)ex.value,
              (long long)dp.value);
        CHECK(feasible(p, dp, 1, 120), "dp schedule infeasible");
    }
    std::printf("dp: 50 instances of 300 passes match\n");

    // 6. Monotone in k, and once k reaches peak concurrency every pass is taken.
    {
        auto p = random_passes(rng, 200, 20000);
        std::vector<int> all;
        int64_t total = 0;
        for (size_t j = 0; j < p.size(); ++j) {
            all.push_back(static_cast<int>(j));
            total += p[j].value;
        }
        int peak = sched::max_overlap(p, all, 60);
        int64_t prev = 0;
        for (int k = 1; k <= peak; ++k) {
            auto ex = sched::schedule_mincost_flow(p, k, 60);
            CHECK(ex.value >= prev, "not monotone at k=%d", k);
            CHECK(feasible(p, ex, k, 60), "k=%d infeasible", k);
            prev = ex.value;
        }
        CHECK(prev == total, "k=peak (%d) should take everything: %lld of %lld", peak, (long long)prev,
              (long long)total);
        auto below = sched::schedule_mincost_flow(p, peak - 1, 60);
        CHECK(below.value < total, "k=peak-1 cannot take everything");
        std::printf("saturation: peak concurrency %d, all %zu passes scheduled at k = %d\n", peak, p.size(), peak);
    }

    if (failures) {
        std::printf("%d check(s) failed\n", failures);
        return 1;
    }
    std::printf("all checks passed\n");
    return 0;
}
