# 016 - Sun-Synchronous Repeat Ground Track Design

**Domain:** Orbital mechanics

## Problem

Earth-observation missions want two things from their orbit at once:

1. **Sun-synchronous**: the orbit plane precesses at the mean Sun's rate
   (0.98565 deg/day) so every pass over a site happens at the same local time.
2. **Repeat ground track**: after exactly N revolutions in D days the track
   lies on top of itself, so the same scene is revisited with the same geometry.

Given (N, D), find the semi-major axis and inclination that satisfy both,
then list every such orbit in an altitude band, work out the order in which
adjacent tracks get filled in, and check the design by flying it through a
numerical J2 propagator for a full cycle.

Dependency-free Rust.

## Approach

First-order J2 secular rates for a circular orbit, `k = J2 (Re/a)^2`:

```
dRAAN/dt = -3/2 n k cos i
dw/dt    =  3/4 n k (5 cos^2 i - 1)
dM/dt    =  n [1 + 3/4 k (3 cos^2 i - 1)]
```

- **SSO condition** fixes `i(a)`: `cos i = -dRAAN_sun / (3/2 n k)`.
- **Repeat condition**: N nodal periods equal D nodal days,
  `(dw/dt + dM/dt) = (N/D) (omega_E - dRAAN/dt)`. With `i = i(a)` this is one
  equation in `a`, monotonic over LEO, solved by bracketed bisection.
- **Sub-cycle**: rev k crosses the node `(k D mod N)` track spacings west of
  the reference. The adjacent track is reached after the smallest d with
  `d N = -/+1 (mod D)`, which gives the days between neighbouring swaths and
  whether the neighbour is east or west. Cross-checked by brute force over
  all coprime N/D with D <= 30.
- **Enumeration**: for each D, N is bracketed by the Keplerian revs/day at the
  band edges; only coprime pairs are solved.
- **Numerical check**: Cartesian two-body + J2 acceleration, RK4 at 5 s.
  Mean elements are not osculating elements (J2 short-period terms put the
  osculating `a` about 9 km above the mean at the node of a near-polar orbit),
  so the initial state is iterated until the orbit-averaged `a` and `i` match
  the design. The propagator then measures the nodal period, the node drift
  rate, and how far the 233rd node lands from the first in Earth-fixed
  longitude.

## Verification (`tests/verify.rs`)

| Check | Reference | Result |
|---|---|---|
| SSO and repeat closure identities | exact | 1e-12 |
| SSO nodal day | 86400 s | < 0.05 s |
| Landsat-8 233/16 | 705 km nominal, 98.2 deg | 699.6 km, 98.19 deg |
| Sentinel-1 175/12 | 693 km, 98.18 deg | 692.8 km, 98.16 deg |
| Sentinel-2 143/10 | 786 km, 98.62 deg | 786.1 km, 98.54 deg |
| J2 lowering vs Keplerian | hand estimate `a k (2/3)(3/4)(2 - 8 cos^2 i)` | 0.05 km |
| Short-period `a` offset at node | `3/2 J2 Re^2/a sin^2 i` = 9.13 km | 0.1 km |
| Landsat adjacent-path interval | 7 days | 7 days |
| Propagated nodal period | design | 1.6e-6 rel |
| Propagated node rate | Sun rate | 7.9e-4 rel (J2^2 terms) |
| Track closure after 16 days | < 1% of 172 km spacing | 0.39 km |

Published altitudes are nominal figures and real orbits are tuned with
higher-order geopotential terms, so a few km and under 0.1 deg is the expected
agreement for first-order theory. The 0.08% node-rate residual in the numerical
run corresponds to about 0.3 deg/year of local-time drift, which is what
operational missions trim out with small inclination biases.

## Run

```
cargo run --release
cargo test --release
```

## Sample output

```
Published missions (N revs / D days):
mission         N   D    a [km]   h [km]   i [deg]   T_N [s]  dx [km]  sub-cycle
Landsat-8     233  16   7077.72   699.58   98.1863   5933.05    172.0  7 d, west
Sentinel-1    175  12   7070.97   692.83   98.1588   5924.57    229.0  5 d, west
Sentinel-2    143  10   7164.26   786.12   98.5446   6041.96    280.2  3 d, west

All repeat SSOs between 600 and 800 km with cycles up to 7 days:
   N   D   h [km]   i [deg]   dx [km]     shift  sub-cycle
 104   7   605.50   97.8093     385.3    2697.4  1 d, west
  ...
  29   2   719.95   98.2696    1381.9    2763.8  1 d, east
 100   7   790.91   98.5647     400.8    2805.3  3 d, west

Numerical J2 propagation of the 233/16 design over one full cycle (RK4, 5 s):
  nodal period : design 5933.0472 s, propagated 5933.0380 s (-1.55e-6 rel)
  node rate    : design 0.98565 deg/day, propagated 0.98487 deg/day (-7.91e-4 rel)
  ground track closure after 233 revs: -0.394 km (track spacing 172.0 km)
```
