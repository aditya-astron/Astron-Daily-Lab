# 031 - Flight Software Schedulability: Response-Time Analysis with the Priority Ceiling Protocol

**Domain:** Embedded firmware

## Problem

A small-satellite flight software build runs seven periodic tasks on one
core under a preemptive fixed-priority RTOS, from a 10 ms attitude loop down
to a 500 ms housekeeping logger. Tasks share three mutex-protected resources
(state database, IMU bus, CAN bus). Total utilisation is 0.89, so the
utilisation-based tests cannot answer whether every deadline is met. We want
the exact worst-case response time of every task including blocking on
shared resources, the WCET margin left in the build, and evidence that the
analysis is right.

| task | T (ms) | D (ms) | C (ms) | critical sections (ms) |
|---|---|---|---|---|
| aocs_ctrl | 10 | 10 | 1.8 | state_db 0.3, imu_bus 0.4 |
| rw_cmd | 20 | 15 | 2.5 | can_bus 0.4 |
| nav_filter | 50 | 50 | 9.0 | imu_bus 0.6, state_db 1.0 |
| fdir | 100 | 80 | 6.0 | state_db 0.5 |
| tm_packetize | 100 | 100 | 15.0 | state_db 1.2, can_bus 0.8 |
| thermal_ctrl | 200 | 200 | 19.0 | state_db 1.5, can_bus 0.5 |
| hk_logger | 500 | 500 | 50.0 | state_db 2.0 |

## Approach

**Priorities.** Deadline-monotonic (shorter D, higher priority), which is
the optimal fixed-priority order when D <= T.

**Locking.** Immediate ceiling priority protocol (ICPP, the `PRIO_PROTECT`
mutex in POSIX): a task that takes a mutex runs at the resource ceiling, the
priority of its highest-priority user. On one core a job then never finds a
lock held once it starts, deadlock is impossible, and a job is blocked at
most once, by one lower-priority critical section:

    B_i = max { len(cs) : cs in a lower-priority task, ceiling(res(cs)) >= prio_i }

**Exact test.** Response-time analysis, iterated to a fixed point from
`R = C_i + B_i`:

    R_i = C_i + B_i + sum_{j in hp(i)} ceil(R_i / T_j) C_j,   schedulable iff R_i <= D_i

**Breakdown scale.** Bisection on a uniform WCET multiplier (critical
sections scaled too) for the largest value that keeps every R_i <= D_i.

**Simulator.** A tick-level (100 us) model of the kernel: per-tick release,
selection of the highest effective priority among ready jobs, and one of
three mutex policies: ICPP, priority inheritance, or a plain mutex with no
priority change. Under ICPP it counts any case where the selected job finds
its lock held, which the theory says never happens.

## Verification (`make test`)

- Textbook three-task set (T,C) = (7,3), (12,3), (20,5): R = 3, 6, 20 by
  hand; adding one tick to the third task makes it unschedulable.
- Hand values for the flight set: ceilings, U, B_aocs = 2.0 ms from the
  logger, R_aocs = 18 + 20 = 38 ticks, R_rw = 25 + 20 + 18 = 63 ticks.
- Critical instant: with synchronous release and no lock held, the first
  simulated job of every task matches RTA with B = 0 exactly.
- Tight blocking: for each task, the worst blocker enters its critical
  section at t = 0 and the task and everything above it are released at
  t = 1. The simulated response equals RTA with B = L - 1 (one tick of the
  section is already done), for all six tasks with nonzero blocking.
- 300 random phasings with random execution times over two hyperperiods:
  no misses, no ICPP lock conflicts, every response below its bound.
- Priority inversion by hand: with the logger holding `state_db` and one job
  of everything else released a tick later, ICPP and inheritance give aocs
  37 ticks; a plain mutex lets rw_cmd (25), nav_filter up to its lock (80)
  and thermal_ctrl up to its lock (170) run first: 25 + 80 + 170 + 19 + 18 =
  312 ticks, which the simulator reproduces.
- The breakdown scale is a boundary of the exact test and lies in (1, 1/U).
- Also built and run clean under AddressSanitizer and UBSan.

## How to run

    make run     # schedulability report
    make test    # checks

Needs a C11 compiler and libm.

## Sample output

```
Tasks 7, utilisation U = 0.8900, hyperperiod 1000 ms
Liu-Layland bound n(2^(1/n)-1) = 0.7286 -> inconclusive
Hyperbolic bound prod(1+Ui) = 2.3000 (<= 2 needed) -> inconclusive
Resource ceilings: state_db=P0 imu_bus=P0 can_bus=P1

task          prio   T ms   D ms   C ms   B ms     R ms blocked by          slack
aocs_ctrl        0   10.0   10.0    1.8    2.0      3.8 hk_logger/state_db    62%
rw_cmd           1   20.0   15.0    2.5    2.0      6.3 hk_logger/state_db    58%
nav_filter       2   50.0   50.0    9.0    2.0     17.1 hk_logger/state_db    66%
fdir             3  100.0   80.0    6.0    2.0     27.4 hk_logger/state_db    66%
tm_packetize     4  100.0  100.0   15.0    2.0     48.5 hk_logger/state_db    52%
thermal_ctrl     5  200.0  200.0   19.0    2.0     88.7 hk_logger/state_db    56%
hk_logger        6  500.0  500.0   50.0    0.0    296.5 -                     41%

Breakdown WCET scale: 1.0889 (every WCET can grow 8.9% before a deadline is missed)

Simulated (200 random phasings x 2 hyperperiods, ICPP): misses 0, lock conflicts 0
  aocs_ctrl     worst observed    3.7 ms   bound    3.8 ms
  rw_cmd        worst observed    5.8 ms   bound    6.3 ms
  nav_filter    worst observed   15.1 ms   bound   17.1 ms
  fdir          worst observed   20.2 ms   bound   27.4 ms
  tm_packetize  worst observed   39.7 ms   bound   48.5 ms
  thermal_ctrl  worst observed   68.3 ms   bound   88.7 ms
  hk_logger     worst observed  154.0 ms   bound  296.5 ms

Inversion case: hk_logger holds state_db at t=0, everything else released once at t=1
  ICPP         aocs_ctrl response    3.7 ms (D = 10.0 ms)
  inheritance  aocs_ctrl response    3.7 ms (D = 10.0 ms)
  plain mutex  aocs_ctrl response   31.2 ms (D = 10.0 ms)  DEADLINE MISSED
```

Each task keeps 41-66 % of its own deadline as slack, yet the build only
tolerates 8.9 % uniform WCET growth: at U = 0.89 the logger's busy period
grows much faster than linearly. The 2 ms logger critical section on
`state_db` sets the blocking term for every other task, so shortening it (or
moving the logger's write out of the shared lock) is the cheapest margin
available. The last block is the Mars Pathfinder failure in miniature: with a
plain mutex, medium-priority work runs while the attitude loop waits on the
logger, and the 10 ms deadline is missed by a factor of three.
