# 035 - Beacon Spiral-Scan Acquisition for a LEO-GEO Optical Relay Link

**Domain:** Laser/optical communication

## Problem

A GEO data-relay terminal must find a LEO user 45 000 km away before any
laser link exists. Its open-loop pointing (ephemeris plus attitude knowledge)
is good to 2.5 mrad per axis, and a 4 W beacon spread over that whole cone
would never close the link. So the relay sweeps a narrow beacon over the
uncertainty region while the user stares with a 13.5 cm acquisition sensor
(NEP 40 fW/rtHz, 10 ms integration, Pfa 1e-9, Pd 0.99). The coarse pointing
assembly slews at most 0.1 rad/s.

What beacon divergence, scan pitch and speed minimise acquisition time while
covering 99.9 % of the error distribution, and does a continuous spiral beat
hexagonal step-and-stare with a 4 ms settle?

## Approach

**Threshold.** For a threshold test on Gaussian noise,
`SNR = Q^-1(Pfa) + Q^-1(1-Pd) = 5.998 + 2.326`, and
`P_th = SNR * NEP * sqrt(1 / (2 tau))` = 2.35 pW.

**Divergence.** A Gaussian far field with 1/e^2 half-angle `w` delivers
`P_rx(theta) = (K / w^2) exp(-2 theta^2 / w^2)`, `K = 2 P A eta / (pi L^2)`.
The detection radius is `theta_d^2 = (w^2/2) ln(K / (w^2 P_th))`;
setting `d theta_d / dw = 0` gives `K / (w^2 P_th) = e`: the best beacon has
exactly 4.34 dB of margin on axis and `theta_d = w / sqrt(2)`.

**Pitch and speed.** A target midway between arms sits `p/2` off the track
and is lit for a chord `2 sqrt(theta_d^2 - p^2/4)`, which must be at least
`v tau`. Maximising the swept-area rate `p v` gives `p = v tau = sqrt(2) theta_d`
when dwell-limited; under a slew limit, fly at `v_max` and open the pitch to
`2 sqrt(theta_d^2 - (v_max tau / 2)^2)`.

**Curvature.** That rule assumes straight arms. Near the centre the track
curls inside the detection circle and the lit arc is shorter: flown at
constant speed, the straight-track pitch leaves holes and has to close to
80 % to cover the disk. Flying only the first turn at 60 % speed recovers
almost all of it. Far out the curvature helps: the outer arm bends toward a
mid-gap target and stretches its dwell by `sqrt(r/(r - p/2)) ~ 1 + p/(4r)`.
The flown pitch is certified by bisection against a dense polar grid of the
whole disk (dwell sampled at `tau/100`) and then backed off 1 %.

**Scan radius and timing.** With per-axis sigma the target radius is
Rayleigh, so `R = sigma sqrt(-2 ln(1 - P_unc))`. The Archimedean spiral
`r = c phi` has exact arc length `(c/2)(phi sqrt(1+phi^2) + asinh phi)`,
inverted by Newton for the beam position at each time step. Covering at rate
`A = p v` makes `T_acq` exponential with rate `A / (2 pi sigma^2)`; the
semi-analytic model averages the exact time-to-radius over the Rayleigh
distribution instead, which includes the slow first turn.

**Step-and-stare.** Spots of radius `theta_d` on a triangular lattice with
spacing `sqrt(3) theta_d` tile the plane with no gaps; each costs
`t_settle + tau` and covers `(3 sqrt3 / 2) theta_d^2`.

**Monte Carlo.** 40 000 true LEO directions drawn from the pointing error;
the user declares detection after the beam has been within `theta_d` for one
full `tau`. Candidate samples per target are restricted to the arms that pass
within `theta_d` (radial band plus an `asin(theta_d / rho)` bearing window).

**Checks** (`tests/verify.rs`):

- `erfc`, `Q^-1(1e-9)`, `Q^-1(0.01)` against tabulated values; threshold and
  received power recomputed by hand.
- Optimal divergence equals a 200 000-point brute-force search; on-axis
  margin is `e`, `theta_d = w/sqrt(2)`, `P_rx(theta_d) = P_th`.
- Closed-form pitch/speed equal a brute-force area-rate maximum in the
  dwell- and slew-limited regimes; worst-case chord equals `v tau`.
- Spiral arc length vs. trapezoid quadrature, Newton inverse round trip,
  beam speed by finite difference in both speed segments.
- Windowed detector identical to a full-track scan on 3000 targets.
- Far-field dwell matches the curved-arc formula to 0.5 % and falls to `tau`.
- Certified spiral covers grids 2x and 4x finer than the certifying grid;
  the constant-speed straight-track spiral does not cover.
- Hex step-and-stare covers the disk; spot count within 6 % of area ratio.
- Monte Carlo `P(acq) >= 0.999` (3 binomial sigmas); mean time within 4 % of
  the semi-analytic spiral model and the exponential step-stare model.

**Limits.** Static target during the scan (LEO-GEO relative angular motion assumed
compensated from ephemeris), no platform jitter or point-ahead, no background
or sun in the sensor FOV, ideal Gaussian far field, and only the first half
of the handshake (the user's return beacon and fine-tracking handover are not
modelled). Step-and-stare settle time is taken as constant per step.

## How to run

    cargo run --release      # scenario (about 15 s, mostly pitch certification)
    cargo test --release     # verification

Rust 2021, no external crates.

## Sample output

```
Beacon acquisition, GEO relay -> LEO user, 45000 km
  required SNR      8.324  (Pfa 1e-9, Pd 0.99)
  threshold power   2.354e-12 W  (-86.28 dBm)

Beacon divergence trade (1/e^2 half-angle w):
     w [urad]   on-axis margin [dB]   theta_d [urad]
       450.0                 12.30            535.5
       675.0                  8.78            678.6
       900.0                  6.28            765.3
      1125.0                  4.34            795.5  <- optimum
      1406.3                  2.40            739.9
      1687.5                  0.82            518.8
      1800.0                  0.26            311.8

Scan design
  pointing sigma    2500 urad/axis, cover 0.999 -> R = 9292 urad
  dwell-limited     v = 0.113 rad/s  p = 1125 urad
  straight track    v = 0.100 rad/s  p = 1237 urad  (slew limited)
  certified         p = 1197 urad (96.8 % of straight), first turn at 60 % speed
  certified, no slow first turn: p = 991 urad (80.1 %)

Acquisition, 40000 Monte Carlo pointing errors
                     area rate      full scan   P(acq)   E[T] model   E[T] MC   T95 MC
  spiral           1.197e-4 rad2/s     2.919 s   0.9999      0.357 s     0.364 s   1.031 s
  step-stare (199)   1.174e-4 rad2/s     2.786 s   0.9999      0.334 s     0.335 s   1.008 s
  spiral model: exact arc-length timing with the slow first turn
  pure area-rate (exponential) model: spiral E[T] 0.328 s, P 0.9999; step-stare P 0.9998
  step-stare area rate equals the spiral's at a settle time of 3.73 ms

Beacon power trade (optimal divergence, same scanner):
   P [W]   w [urad]   theta_d [urad]   regime   E[T] [s]   full scan [s]
     0.5        398              281   dwell       2.605          19.419
     1.0        563              398   dwell       1.317          10.048
     2.0        795              563   dwell       0.669           5.254
     4.0       1125              795   slew        0.357           2.919
     8.0       1591             1125   slew        0.242           2.068
    16.0       2250             1591   slew        0.187           1.664
```

The 4 ms settle sits just above the 3.7 ms break-even, so on area rate alone
the spiral wins by about 2 %, yet step-and-stare still acquires faster on
average: it spends no extra time on a slow first turn, and that turn sits
exactly where the Rayleigh distribution puts most targets. Above 4 W the scan
is slew-limited and more beacon power buys less and less.
