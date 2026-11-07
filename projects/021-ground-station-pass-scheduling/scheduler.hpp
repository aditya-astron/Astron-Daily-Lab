// Astron Built this AP
// Project 021 | Algorithmic problem solving | Ground-station antenna scheduling by min-cost flow
// Author: Aditya (aditya-astron)
// Date: 2026-11-07 IST
//
// Pass windows, a pass-geometry model for circular LEO orbits, and three
// schedulers for k identical antennas at one ground station:
//   - exact:   min-cost flow on the time-line graph (Arkin & Silverberg 1987)
//   - dp:      weighted interval scheduling, exact for k = 1 only
//   - greedy:  highest value first, the usual operations heuristic
#pragma once

#include <cstdint>
#include <vector>

namespace sched {

struct Pass {
    int sat = 0;          // satellite index
    int64_t aos = 0;      // acquisition of signal [s]
    int64_t los = 0;      // loss of signal [s]
    int64_t value = 0;    // downlinkable volume [kbit], integer keeps the LP exact
};

struct Schedule {
    int64_t value = 0;               // total scheduled volume [kbit]
    std::vector<int> antenna;        // per pass: antenna index, or -1 if not tracked
};

// Two passes conflict on one antenna if the second AOS falls before the first
// LOS plus the turnaround (slew back to the horizon, re-point, re-lock).
inline int64_t busy_end(const Pass& p, int64_t turnaround) { return p.los + turnaround; }

// Largest number of passes simultaneously occupying antennas. For interval
// graphs this equals the chromatic number, so a subset is feasible on k
// antennas iff this is <= k.
int max_overlap(const std::vector<Pass>& passes, const std::vector<int>& subset, int64_t turnaround);

// Interval-graph colouring: assigns each selected pass to an antenna using the
// earliest-free rule. Uses exactly max_overlap antennas.
std::vector<int> assign_antennas(const std::vector<Pass>& passes, const std::vector<int>& subset,
                                 int64_t turnaround);

Schedule schedule_mincost_flow(const std::vector<Pass>& passes, int k, int64_t turnaround);
Schedule schedule_dp_single(const std::vector<Pass>& passes, int64_t turnaround);
Schedule schedule_greedy(const std::vector<Pass>& passes, int k, int64_t turnaround);

// Circular-orbit pass geometry over a spherical, non-rotating Earth.
// min_elev_rad: station mask; offset_rad: ground-track central angle from the
// station at closest approach. Returns pass duration in seconds (0 if the
// pass never clears the mask) and the peak elevation.
struct PassGeometry {
    double duration_s;
    double max_elev_rad;
};
PassGeometry pass_geometry(double alt_km, double min_elev_rad, double offset_rad);
double earth_central_angle_limit(double alt_km, double min_elev_rad);

}  // namespace sched
