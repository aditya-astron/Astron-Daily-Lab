// Astron Built this AP
// Project 021 | Algorithmic problem solving | Ground-station antenna scheduling by min-cost flow
// Author: Aditya (aditya-astron)
// Date: 2026-11-07 IST

#include "scheduler.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <queue>
#include <stdexcept>

namespace sched {

namespace {

constexpr double kMuEarth = 398600.4418;  // km^3/s^2
constexpr double kReEarth = 6378.137;     // km
constexpr int64_t kInf = std::numeric_limits<int64_t>::max() / 4;

// Plain successive-shortest-path min-cost flow. Costs can be negative on the
// forward graph, but every edge goes from an earlier to a later time node, so
// the graph is a DAG and one topological relaxation gives valid potentials.
// After that, reduced costs are non-negative and Dijkstra is safe.
class MinCostFlow {
public:
    explicit MinCostFlow(int n) : g_(n) {}

    int add_edge(int u, int v, int64_t cap, int64_t cost) {
        g_[u].push_back({v, static_cast<int>(g_[v].size()), cap, cost});
        g_[v].push_back({u, static_cast<int>(g_[u].size()) - 1, 0, -cost});
        return static_cast<int>(g_[u].size()) - 1;
    }

    // Requires node indices to be a topological order of the forward graph.
    int64_t run(int s, int t, int64_t want, int64_t& flow_out) {
        const int n = static_cast<int>(g_.size());
        std::vector<int64_t> pot(n, kInf);
        pot[s] = 0;
        for (int u = 0; u < n; ++u) {
            if (pot[u] == kInf) continue;
            for (const Edge& e : g_[u])
                if (e.cap > 0 && pot[u] + e.cost < pot[e.to]) pot[e.to] = pot[u] + e.cost;
        }

        int64_t flow = 0, cost = 0;
        std::vector<int64_t> dist(n);
        std::vector<int> prev_node(n), prev_edge(n);
        using QE = std::pair<int64_t, int>;
        while (flow < want) {
            std::fill(dist.begin(), dist.end(), kInf);
            dist[s] = 0;
            std::priority_queue<QE, std::vector<QE>, std::greater<QE>> pq;
            pq.push({0, s});
            while (!pq.empty()) {
                auto [d, u] = pq.top();
                pq.pop();
                if (d > dist[u]) continue;
                for (int i = 0; i < static_cast<int>(g_[u].size()); ++i) {
                    const Edge& e = g_[u][i];
                    if (e.cap <= 0 || pot[e.to] == kInf) continue;
                    int64_t nd = d + e.cost + pot[u] - pot[e.to];
                    if (nd < dist[e.to]) {
                        dist[e.to] = nd;
                        prev_node[e.to] = u;
                        prev_edge[e.to] = i;
                        pq.push({nd, e.to});
                    }
                }
            }
            if (dist[t] == kInf) break;
            for (int v = 0; v < n; ++v)
                if (dist[v] < kInf) pot[v] += dist[v];

            int64_t push = want - flow;
            for (int v = t; v != s; v = prev_node[v])
                push = std::min(push, g_[prev_node[v]][prev_edge[v]].cap);
            for (int v = t; v != s; v = prev_node[v]) {
                Edge& e = g_[prev_node[v]][prev_edge[v]];
                e.cap -= push;
                g_[v][e.rev].cap += push;
                cost += push * e.cost;
            }
            flow += push;
        }
        flow_out = flow;
        return cost;
    }

    int64_t residual(int u, int idx) const { return g_[u][idx].cap; }

private:
    struct Edge {
        int to;
        int rev;
        int64_t cap;
        int64_t cost;
    };
    std::vector<std::vector<Edge>> g_;
};

std::vector<int> all_indices(size_t n) {
    std::vector<int> v(n);
    std::iota(v.begin(), v.end(), 0);
    return v;
}

Schedule finish(const std::vector<Pass>& passes, const std::vector<int>& chosen, int64_t turnaround) {
    Schedule s;
    s.antenna = assign_antennas(passes, chosen, turnaround);
    for (int j : chosen) s.value += passes[j].value;
    return s;
}

}  // namespace

int max_overlap(const std::vector<Pass>& passes, const std::vector<int>& subset, int64_t turnaround) {
    // Sweep over +1 at AOS, -1 at busy end. Ends sort before starts at the same
    // instant because intervals are half-open.
    std::vector<std::pair<int64_t, int>> ev;
    ev.reserve(2 * subset.size());
    for (int j : subset) {
        ev.push_back({passes[j].aos, +1});
        ev.push_back({busy_end(passes[j], turnaround), -1});
    }
    std::sort(ev.begin(), ev.end());
    int cur = 0, best = 0;
    for (auto& e : ev) best = std::max(best, cur += e.second);
    return best;
}

std::vector<int> assign_antennas(const std::vector<Pass>& passes, const std::vector<int>& subset,
                                 int64_t turnaround) {
    std::vector<int> order = subset;
    std::sort(order.begin(), order.end(), [&](int a, int b) { return passes[a].aos < passes[b].aos; });

    std::vector<int> antenna(passes.size(), -1);
    using Slot = std::pair<int64_t, int>;  // (free from, antenna id)
    std::priority_queue<Slot, std::vector<Slot>, std::greater<Slot>> busy;
    std::vector<int> idle;                 // freed antenna ids, lowest reused first
    int next_id = 0;
    for (int j : order) {
        while (!busy.empty() && busy.top().first <= passes[j].aos) {
            idle.push_back(busy.top().second);
            busy.pop();
        }
        int id;
        if (!idle.empty()) {
            auto it = std::min_element(idle.begin(), idle.end());
            id = *it;
            idle.erase(it);
        } else {
            id = next_id++;
        }
        antenna[j] = id;
        busy.push({busy_end(passes[j], turnaround), id});
    }
    return antenna;
}

Schedule schedule_mincost_flow(const std::vector<Pass>& passes, int k, int64_t turnaround) {
    if (k < 1) throw std::invalid_argument("need at least one antenna");
    if (passes.empty()) return Schedule{0, {}};

    // Nodes are the distinct event times. A chain t_i -> t_{i+1} with capacity
    // k carries idle antennas; each pass is a capacity-1 shortcut from its AOS
    // to its busy end with cost -value. Any integral k-unit flow decomposes
    // into k antenna timelines of non-overlapping passes, and conversely, so
    // the min-cost flow is the max-value schedule. The constraint matrix is
    // totally unimodular, which is why integer optimality comes for free.
    std::vector<int64_t> times;
    for (const Pass& p : passes) {
        times.push_back(p.aos);
        times.push_back(busy_end(p, turnaround));
    }
    std::sort(times.begin(), times.end());
    times.erase(std::unique(times.begin(), times.end()), times.end());
    auto node = [&](int64_t t) {
        return static_cast<int>(std::lower_bound(times.begin(), times.end(), t) - times.begin());
    };

    const int n = static_cast<int>(times.size());
    MinCostFlow mcf(n);
    for (int i = 0; i + 1 < n; ++i) mcf.add_edge(i, i + 1, k, 0);
    std::vector<std::pair<int, int>> pass_edge(passes.size());
    for (size_t j = 0; j < passes.size(); ++j) {
        int u = node(passes[j].aos);
        int v = node(busy_end(passes[j], turnaround));
        pass_edge[j] = {u, mcf.add_edge(u, v, 1, -passes[j].value)};
    }

    int64_t flow = 0;
    int64_t cost = mcf.run(0, n - 1, k, flow);

    std::vector<int> chosen;
    for (size_t j = 0; j < passes.size(); ++j)
        if (mcf.residual(pass_edge[j].first, pass_edge[j].second) == 0) chosen.push_back(static_cast<int>(j));

    Schedule s = finish(passes, chosen, turnaround);
    if (s.value != -cost) throw std::logic_error("flow cost and recovered schedule disagree");
    return s;
}

Schedule schedule_dp_single(const std::vector<Pass>& passes, int64_t turnaround) {
    // Classic weighted interval scheduling: sort by busy end, and for each pass
    // either skip it or take it plus the best schedule ending by its AOS.
    const int n = static_cast<int>(passes.size());
    std::vector<int> order = all_indices(passes.size());
    std::sort(order.begin(), order.end(),
              [&](int a, int b) { return busy_end(passes[a], turnaround) < busy_end(passes[b], turnaround); });
    std::vector<int64_t> ends(n);
    for (int i = 0; i < n; ++i) ends[i] = busy_end(passes[order[i]], turnaround);

    std::vector<int64_t> best(n + 1, 0);
    std::vector<int> pred(n);
    for (int i = 0; i < n; ++i) {
        const Pass& p = passes[order[i]];
        pred[i] = static_cast<int>(std::upper_bound(ends.begin(), ends.begin() + i, p.aos) - ends.begin());
        best[i + 1] = std::max(best[i], best[pred[i]] + p.value);
    }
    std::vector<int> chosen;
    for (int i = n; i > 0;) {
        if (best[i] == best[i - 1]) {
            --i;  // pass i-1 not needed for this optimum
        } else {
            chosen.push_back(order[i - 1]);
            i = pred[i - 1];
        }
    }
    return finish(passes, chosen, turnaround);
}

Schedule schedule_greedy(const std::vector<Pass>& passes, int k, int64_t turnaround) {
    std::vector<int> order = all_indices(passes.size());
    std::stable_sort(order.begin(), order.end(),
                     [&](int a, int b) { return passes[a].value > passes[b].value; });
    std::vector<int> chosen;
    for (int j : order) {
        chosen.push_back(j);
        if (max_overlap(passes, chosen, turnaround) > k) chosen.pop_back();
    }
    return finish(passes, chosen, turnaround);
}

double earth_central_angle_limit(double alt_km, double min_elev_rad) {
    // Largest station-to-subsatellite central angle at which the satellite is
    // still above the elevation mask: lambda = acos(Re cos(eps) / r) - eps.
    const double r = kReEarth + alt_km;
    return std::acos(kReEarth * std::cos(min_elev_rad) / r) - min_elev_rad;
}

PassGeometry pass_geometry(double alt_km, double min_elev_rad, double offset_rad) {
    const double r = kReEarth + alt_km;
    const double lam_max = earth_central_angle_limit(alt_km, min_elev_rad);
    if (std::fabs(offset_rad) >= lam_max) return {0.0, 0.0};

    // Right spherical triangle: station, closest-approach point, and the point
    // where the track crosses the mask circle. cos(lam_max) = cos(off) cos(half).
    const double half_arc = std::acos(std::cos(lam_max) / std::cos(offset_rad));
    const double n = std::sqrt(kMuEarth / (r * r * r));
    const double duration = 2.0 * half_arc / n;

    // Elevation at central angle lambda: tan(el) = (cos(lam) - Re/r) / sin(lam).
    const double lam = std::fabs(offset_rad);
    const double el = lam < 1e-12 ? M_PI / 2
                                  : std::atan2(std::cos(lam) - kReEarth / r, std::sin(lam));
    return {duration, el};
}

}  // namespace sched
