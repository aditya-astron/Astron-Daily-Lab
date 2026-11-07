# 019 - Mars X-band Downlink Over a Synodic Cycle

**Domain:** RF link budgets

## Problem

A Mars orbiter with a 100 W X-band TWTA and a 3 m high-gain antenna downlinks
to a DSN 34 m beam-waveguide antenna at Goldstone. Earth-Mars range swings
from 0.68 AU to 2.4 AU over the next synodic period, so the supportable data
rate changes by more than a factor of ten, and around solar conjunction the
link is lost entirely. How much data can one 8 h track per day return across
the full cycle starting 2026-11-07, when is the outage, and how much does
stepping the data rate through the pass (variable data rate) buy over a single
fixed rate?

Python with numpy/scipy.

## Approach

**Geometry** (`ephem.py`). Earth and Mars from the JPL approximate Keplerian
elements with secular rates (Standish, valid 1800-2050). Kepler's equation by
Newton. From the heliocentric vectors: range, one-way light time, Sun-Earth-Probe
angle, and geocentric declination (ecliptic to equatorial at the J2000 obliquity).

**Link budget** (`link.py`):

```
Pr/N0 = Pt + Gt - Ltx - Lpt - Lfs - Latm(el) - Lpol + Gr - Lpt_gs - 10 log10(k Tsys(el))
Lfs   = 20 log10(4 pi d / lambda)
Gt    = eta (pi D / lambda)^2
Latm  = A_zenith / sin(el)
Tsys  = T_rx + T_phys (1 - 1/Latm) + T_cmb / Latm
Rb    = Pr/N0 - (Eb/N0)_req - margin          [dB-bps], floored onto a 2^k rate ladder
```

**Coding threshold.** The required Eb/N0 is anchored to a computed limit rather
than a quoted number: the binary-input AWGN capacity

```
C(Es/N0) = 1 - E[ log2(1 + exp(-2y/s^2)) ],  y ~ N(1, s^2),  s^2 = 1 / (2 Es/N0)
```

is evaluated by Gauss-Hermite quadrature and solved for `C = R` at `R = 1/2`,
giving 0.187 dB. The design threshold adds a 1.0 dB gap for a practical
finite-length turbo/LDPC code, then a 3 dB link margin on top.

**Passes.** For each day the track is centred on transit, sampled every
2 minutes with `sin el = sin(lat) sin(dec) + cos(lat) cos(dec) cos(H)`, and
counted only above 10 deg. VDR uses the best ladder rate at each sample's
elevation; fixed rate uses the rate that closes at the lowest elevation of the
track. Days with SEP below 3 deg are treated as conjunction outage.

## Verification (`test_mars_link.py`, 17 checks)

| Check | Reference | Result |
|---|---|---|
| Julian date | J2000.0 and Meeus ex. 7.a | exact |
| Kepler residual, e up to 0.97 | 0 | < 1e-12 rad |
| Earth perihelion / aphelion | 0.98329 / 1.01671 AU, early January | within 2e-4 AU |
| Mars closest approach 2025-01-12 | 96.08 million km | 96.09, same day |
| Mars opposition 2025-01-16 | published date | same day |
| Mars solar conjunction 2026-01-09 | published date | same day, SEP 0.94 deg |
| Mars closest approach 2027-02-20 | 101.4 million km | 101.40, same day |
| Uncoded BPSK at BER 1e-5 / 1e-6 | 9.59 / 10.53 dB | 9.588 / 10.530 dB |
| BPSK Monte Carlo, 2M bits at 6 dB | theory | within 0.05 dB |
| BI-AWGN limit, R = 1/2, 1/3, 1/6 | 0.187, -0.495, -1.08 dB | 0.187, -0.495, -1.073 dB |
| Shannon limit R -> 0 | -1.59 dB | exact |
| Free-space loss, 1 AU at 8.42 GHz | hand calc 274.45 dB | 1e-9 dB |
| 34 m gain at 74 % aperture efficiency | 68.24 dBi | 0.02 dB |
| Rate vs range | inverse square | exact |
| Fixed-rate pass volume | hand calc | exact |

## Run

```
python3 mars_link.py                    # cycle from 2026-11-07
python3 mars_link.py 2028-01-01 120     # any start date and span
python3 -m pytest -q                    # or: python3 test_mars_link.py
```

## Sample output

```
Link budget at farthest approach, 2028-05-12, OWLT 20.1 min, peak el 69 deg:
  EIRP                                       64.56 dBW
  Free-space loss (2.412 AU)               -282.10 dB
  Tsys (22.4 K)                              13.50 dBK
  Pr/N0                                      65.55 dB-Hz
  Max information rate                       61.36 dB-bps
  -> rate ladder                              1.024 Mbps (unquantised 1.368 Mbps)

 month     range AU     SEP deg   peak rate     Gbit/day VDR   Gbit/day fixed
 2027-02  0.678-0.718   154.1    8.192 Mbps        235.93          235.93
 2027-08  1.775-1.941    51.3    2.048 Mbps         49.57           33.30
 2027-12  2.243-2.295    18.1    1.024 Mbps         27.75           14.06
 2028-03  2.370-2.396     0.8    1.024 Mbps          6.66            6.66

Solar conjunction 2028-03-21, min SEP 0.81 deg
Outage (SEP < 3 deg): 2028-03-08 to 2028-04-03, 27 days
Cycle volume, one 8 h track/day: VDR 7.40 TB, fixed rate 6.96 TB, VDR gain 6.3 %
Best/worst daily volume ratio: 8.6x  (range^2 ratio 12.7x)
```

## Notes

- Near opposition the link closes at about 17 Mbps but the rate table tops out
  at 8.192 Mbps, so volume flattens there. That cap, not physics, is why the
  best/worst volume ratio (8.6x) is below the range-squared ratio (12.7x).
- VDR pays most when Mars sits far south of Goldstone's latitude (December
  2027: 2x), because the fixed rate is then pinned by a low-elevation track
  edge. With a power-of-two ladder the cycle-wide gain is modest; a finer
  ladder recovers more.
- The rates land in the same range as flown Mars orbiters with similar
  hardware (hundreds of kbps to a few Mbps into a 34 m antenna).
- Not modelled: solar noise rise and scintillation between 3 and ~10 deg SEP,
  rain fade beyond the clear-sky zenith value, and residual-carrier power
  split for ranging.
