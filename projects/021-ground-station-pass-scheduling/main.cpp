// Astron Built this AP
// Project 021 | Algorithmic problem solving | Ground-station antenna scheduling by min-cost flow
// Author: Aditya (aditya-astron)
// Date: 2026-11-07 IST
//
// One day of LEO passes over a single mid-latitude ground station with 1 to 3
// identical antennas. Pass durations come from circular-orbit geometry over a
// 10 deg mask; each satellite has its own downlink rate, so a pass is worth
// rate x duration. Compares the exact min-cost-flow schedule with the
// highest-value-first greedy rule that is common in operations tooling.

#include <algorithm>
#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <vector>

#include "scheduler.hpp"

namespace {

// Small deterministic generator so the demo is reproducible across compilers.
struct Lcg {
    uint64_t s;
    double uniform() {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        return static_cast<double>(s >> 11) * (1.0 / 9007199254740992.0);
    }
};

struct Sat {
    double alt_km;
    double rate_kbps;
};

std::vector<sched::Pass> one_day_of_passes(const std::vector<Sat>& sats, double min_elev) {
    constexpr double kDay = 86400.0;
    constexpr double kSidereal = 86164.1;
    constexpr double kStationLatCos = 0.766;  // ~40 deg latitude
    Lcg rng{20261107};
    std::vector<sched::Pass> out;
    for (int s = 0; s < static_cast<int>(sats.size()); ++s) {
        const double r = 6378.137 + sats[s].alt_km;
        const double period = 2 * M_PI * std::sqrt(r * r * r / 398600.4418);
        // Ground-track shift between consecutive revs, measured at the station.
        const double shift = 2 * M_PI * period / kSidereal * kStationLatCos;
        const double lam_max = sched::earth_central_angle_limit(sats[s].alt_km, min_elev);
        // Two visibility windows per day (ascending and descending node sides),
        // roughly 12 h apart, each a run of consecutive revolutions whose
        // closest-approach offset sweeps across the station.
        for (int side = 0; side < 2; ++side) {
            double t_mid = (side * 0.5 + 0.45 * rng.uniform()) * kDay;
            double offset = -lam_max + 0.6 * shift * rng.uniform();
            for (; offset < lam_max; offset += shift, t_mid += period) {
                sched::PassGeometry g = sched::pass_geometry(sats[s].alt_km, min_elev, offset);
                if (g.duration_s < 60.0) continue;
                if (t_mid + 0.5 * g.duration_s > kDay) break;  // keep the plan inside one day
                sched::Pass p;
                p.sat = s;
                p.aos = std::llround(t_mid - 0.5 * g.duration_s);
                p.los = std::llround(t_mid + 0.5 * g.duration_s);
                p.value = std::llround(sats[s].rate_kbps * g.duration_s);
                out.push_back(p);
            }
        }
    }
    return out;
}

void hhmm(int64_t t, char* buf) {
    std::snprintf(buf, 8, "%02d:%02d", static_cast<int>(t / 3600) % 24, static_cast<int>(t / 60) % 60);
}

}  // namespace

int main() {
    const double min_elev = 10.0 * M_PI / 180.0;
    const int64_t turnaround = 120;  // s: slew to next AOS azimuth, re-lock, reconfigure
    // 32-satellite mixed fleet: altitudes 480-720 km, downlink 50-600 Mbit/s.
    std::vector<Sat> sats;
    Lcg fleet_rng{42};
    for (int s = 0; s < 32; ++s)
        sats.push_back({480.0 + 240.0 * fleet_rng.uniform(), 1e3 * (50.0 + 550.0 * fleet_rng.uniform())});

    std::vector<sched::Pass> passes = one_day_of_passes(sats, min_elev);
    int64_t demand = 0;
    std::vector<int> every;
    for (size_t j = 0; j < passes.size(); ++j) {
        demand += passes[j].value;
        every.push_back(static_cast<int>(j));
    }

    std::printf("Ground station: 1 site, 10 deg mask, %" PRId64 " s turnaround\n", turnaround);
    std::printf("Fleet: %zu satellites, %zu passes in 24 h, peak concurrency %d\n", sats.size(),
                passes.size(), sched::max_overlap(passes, every, turnaround));
    std::printf("Total offered volume: %.1f Gbit\n\n", demand / 1e6);

    std::printf(" k |  exact (flow)  |     greedy     |  gap   | passes exact/greedy | captured\n");
    std::printf("---+----------------+----------------+--------+---------------------+---------\n");
    for (int k = 1; k <= 3; ++k) {
        sched::Schedule ex = sched::schedule_mincost_flow(passes, k, turnaround);
        sched::Schedule gr = sched::schedule_greedy(passes, k, turnaround);
        int n_ex = 0, n_gr = 0;
        for (int a : ex.antenna) n_ex += a >= 0;
        for (int a : gr.antenna) n_gr += a >= 0;
        std::printf(" %d | %9.1f Gbit | %9.1f Gbit | %5.2f%% | %9d / %-7d | %6.1f%%\n", k, ex.value / 1e6,
                    gr.value / 1e6, 100.0 * (ex.value - gr.value) / ex.value, n_ex, n_gr,
                    100.0 * ex.value / demand);
        if (k == 1) {
            sched::Schedule dp = sched::schedule_dp_single(passes, turnaround);
            std::printf("   | DP cross-check for k = 1: %.1f Gbit (%s)\n", dp.value / 1e6,
                        dp.value == ex.value ? "match" : "MISMATCH");
        }
    }

    // Antenna timeline for the two-antenna case, first 12 tracked passes.
    sched::Schedule two = sched::schedule_mincost_flow(passes, 2, turnaround);
    std::vector<int> order;
    for (size_t j = 0; j < passes.size(); ++j)
        if (two.antenna[j] >= 0) order.push_back(static_cast<int>(j));
    std::sort(order.begin(), order.end(), [&](int a, int b) { return passes[a].aos < passes[b].aos; });
    std::printf("\nTwo-antenna plan (first 12 tracks):\n  ant  sat   AOS    LOS   dur[s]  volume[Gbit]\n");
    for (size_t i = 0; i < order.size() && i < 12; ++i) {
        const sched::Pass& p = passes[order[i]];
        char a[8], b[8];
        hhmm(p.aos, a);
        hhmm(p.los, b);
        std::printf("  %3d  %3d  %s  %s  %6" PRId64 "  %10.2f\n", two.antenna[order[i]], p.sat, a, b,
                    p.los - p.aos, p.value / 1e6);
    }
    return 0;
}
