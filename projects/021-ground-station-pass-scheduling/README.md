# 021 - Ground-Station Antenna Scheduling by Min-Cost Flow

**Domain:** Algorithmic problem solving

## Problem

One ground station with k identical antennas serves a 32-satellite LEO fleet.
Over a day the station sees about 120 passes. Each pass has an AOS, an LOS and
a value (downlink rate x duration), and an antenna needs a fixed turnaround
after LOS to slew and re-lock before it can take the next AOS. Which passes
should be tracked, and on which antenna, to maximise the data brought down?

Operations tools often answer this with "highest value first". That rule is
not optimal even for one antenna, and the question is how much it leaves on
the table and what the exact answer costs to compute.

## Approach

**Pass geometry.** Circular orbit, spherical non-rotating Earth, elevation
mask eps. The largest station-to-subsatellite central angle still above the
mask, and the pass length for a track whose closest approach is offset
`lambda_min` from the station:

```
lambda_max = acos(Re cos(eps) / r) - eps
half_arc   = acos(cos(lambda_max) / cos(lambda_min))      (right spherical triangle)
T_pass     = 2 half_arc / n,   n = sqrt(mu / r^3)
```

Successive revolutions shift the track by roughly `2 pi P / T_sidereal x cos(lat)`,
which produces the familiar run of short-long-short passes in each window.

**Exact schedule (min-cost flow).** Make a node for each distinct event time
(AOS and LOS + turnaround). Chain consecutive nodes with capacity k and cost 0;
these edges carry idle antennas. Add one edge per pass from its AOS node to its
busy-end node with capacity 1 and cost -value. Push k units from the first node
to the last. An integral k-unit flow is exactly k antenna timelines of
non-overlapping passes, and the constraint matrix is totally unimodular, so
the min-cost flow is the max-value schedule (Arkin & Silverberg, 1987).

The graph is a DAG in time order, so one topological relaxation gives valid
potentials despite the negative costs, and each augmentation after that is a
Dijkstra on reduced costs: O(k E log V) overall.

Selected passes are then put on antennas with the earliest-free interval
colouring rule. Interval graphs are perfect, so the number of antennas used
equals the peak concurrency of the selected passes, which the flow keeps <= k.

**Baselines.** Greedy by value (add a pass if peak concurrency stays <= k),
and for k = 1 the classic weighted interval scheduling DP.

## Verification (`make test`)

- Overhead pass at 500 km: 693.3 s with a 0 deg mask and 443.0 s with a
  10 deg mask, both computed by hand from the formulas above; a grazing pass
  has zero length and peaks exactly at the mask.
- Hand-worked case where greedy takes a 10 and the optimum is 7 + 6 = 13,
  plus the same case with turnaround and with two antennas.
- 600 random instances (6-14 passes, k = 1..3, turnaround 0/60/120 s) against
  exhaustive subset search. All match; greedy is strictly worse in 79.
- 50 instances of 300 passes against the k = 1 DP.
- Every schedule is checked for double-booking within the turnaround, antenna
  index < k, and value consistency.
- Optimal value is monotone in k, takes every pass once k equals the peak
  concurrency, and cannot at k - 1.

## Run

```
make
./pass_scheduler
make test
```

Needs a C++17 compiler; no dependencies.

## Sample output

```
Ground station: 1 site, 10 deg mask, 120 s turnaround
Fleet: 32 satellites, 122 passes in 24 h, peak concurrency 3
Total offered volume: 14226.0 Gbit

 k |  exact (flow)  |     greedy     |  gap   | passes exact/greedy | captured
---+----------------+----------------+--------+---------------------+---------
 1 |    9109.9 Gbit |    8971.4 Gbit |  1.52% |        72 / 71      |   64.0%
   | DP cross-check for k = 1: 9109.9 Gbit (match)
 2 |   13160.6 Gbit |   13160.6 Gbit |  0.00% |       111 / 111     |   92.5%
 3 |   14226.0 Gbit |   14226.0 Gbit |  0.00% |       122 / 122     |  100.0%

Two-antenna plan (first 12 tracks):
  ant  sat   AOS    LOS   dur[s]  volume[Gbit]
    0    1  00:18  00:22     242       95.95
    1   27  00:23  00:28     334       96.53
    0   23  00:38  00:40     153        9.48
    0    4  01:01  01:09     455       28.98
    ...
```

With one antenna the greedy rule loses 138 Gbit a day, about two full passes
of the faster satellites. A second antenna recovers 92.5% of the offered
volume and a third, matching the peak concurrency, takes everything. Where
contention is light the greedy rule happens to be optimal, which is why it
survives in practice; the flow solution costs milliseconds and removes the
guesswork.

## Limitations

Single station, fixed pass values, no per-satellite minimum contact
requirement and no rotating-Earth pass geometry. Multi-station scheduling with
per-satellite constraints is no longer a pure flow problem and needs an ILP.
