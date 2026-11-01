# 003 - Lambert Targeting for a LEO Rendezvous Transfer

**Domain:** Orbital mechanics

## Problem

A chaser in a 300 km circular parking orbit has to meet a target in an 800 km
circular orbit (same plane). At t = 0 the target leads by 40 deg. When should
the chaser burn, how long should the transfer take, and what does it cost?
How much more does it cost to get there faster?

The tool for "go from r1 to r2 in time t" is Lambert's problem. This project
writes a solver, checks it against textbook cases and an independent
integrator, then uses it to search the (coast, time-of-flight) space. The
search has no knowledge of Hohmann transfers, so recovering the Hohmann
answer is the end-to-end check.

## Approach

Universal-variable Lambert (Bate-Mueller-White, Curtis Alg. 5.2), single rev:

```
A    = sin(dth) sqrt(r1 r2 / (1 - cos dth))
y(z) = r1 + r2 + A (z S(z) - 1) / sqrt(C(z))
F(z) = (y/C)^1.5 S + A sqrt(y) - sqrt(mu) t = 0
f = 1 - y/r1,  g = A sqrt(y/mu),  gdot = 1 - y/r2
v1 = (r2 - f r1)/g,   v2 = (gdot r2 - r1)/g
```

- `F(z)` is monotonic on the single-rev branch, so it is solved with a
  bracketed Brent root (bracket: just inside `y = 0` on the left, `z -> 4 pi^2`
  on the right) instead of a Newton loop that can wander into `y < 0`.
- Stumpff `C`, `S` use a 6-term Taylor series for `|z| < 0.1`; the closed
  forms lose half their digits near `z = 0`.
- Transfers within 2e-3 rad of 180 deg are refused. There `g -> 0` and the
  plane is undefined; in testing the propagated miss grew from 0.1 mm to
  10 m as the angle moved from 3e-3 to 3e-5 rad off 180 deg.

Rendezvous search: coarse grid over one synodic period of coast and
0.3 - 1.6 Hohmann times of flight, then Nelder-Mead polish (derivative-free,
since the optimum sits on the edge of the refused 180 deg band). Burns are
`|v1 - v_chaser|` and `|v_target - v2|`.

## Verification

- Curtis Example 5.2 and Vallado Example 7-5 velocities.
- Hyperbolic, elliptic and long-way transfers propagated with DOP853
  (rtol 1e-12): endpoint miss < 1 m, velocity error < 1e-7 km/s.
- Energy and angular momentum equal at both ends.
- Hohmann 300 -> 800 km by hand: 138.16 + 135.68 m/s, 47.82 min.
- Lambert search optimum lands within 0.5 m/s, 1 min of the Hohmann answer
  and never below it.

## Run

```
pip install numpy scipy pytest
python rendezvous.py
python -m pytest -q
```

## Sample output

```
Hohmann reference
  dv1  138.16 m/s  dv2  135.68 m/s  total  273.84 m/s
  transfer  47.82 min, coast before burn   75.15 min

Lambert optimum (grid + Nelder-Mead)
  coast   75.14 min  tof  47.86 min
  dv1  138.16 m/s  dv2  135.68 m/s  total  273.84 m/s
  excess over Hohmann 0.001 m/s
  propagated miss distance 3.790e-03 m

Fast-intercept penalty (best coast for each tof)
   tof [min]   best dv [m/s]
      19.01         807.7
      28.85         479.6
      38.18         329.5
      48.02         277.1
      57.35         326.5
```

Cutting the transfer to 40% of Hohmann time costs about 3x the delta-v.
