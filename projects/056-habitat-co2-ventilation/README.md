# 056 - Multi-Zone CO2 Ventilation Network for a Sealed Habitat

**Domain:** Analog habitat systems / ECLSS

## Problem

A four-crew sealed analog habitat has five zones: Hab (120 m^3, holds the
CO2 scrubber), a central Node (30 m^3), Lab (80 m^3), Greenhouse (60 m^3)
and Crew Quarters (36 m^3). Every module trades air with the Node through a
balanced inter-module ventilation (IMV) fan pair. All mixing is mechanical
(doors closed, as in a flight vehicle where buoyancy does nothing), so the
scrubber only sees CO2 that the fans carry to it.

A single well-mixed number says the habitat sits at 1.4 mmHg. The crew
sleeps in the quarters, though, two fan hops away from the scrubber. The
questions:

1. How big must the quarters fan pair be to keep the sleeping crew under
   3.0 mmHg ppCO2 (about 3900 ppm), and does sizing on the steady state
   get it right?
2. If that fan pair fails at lights-out, how long until the limit?
3. If the scrubber fails, which zone hits the limit first, and when?

## Approach

Each zone is well mixed. With `c` the CO2 density in kg/m^3:

```
V_i dc_i/dt = sum_j Q_ij (c_j - c_i) - eta Qs c_i [scrubber zone] + S_i(t)
dc/dt = A c + s,   A = V^-1 L
```

Balanced pairs make `L` symmetric, so `A` is similar to
`V^-1/2 L V^-1/2`: every zone mode is real and decaying. The eigenvalues come
from a Jacobi solve. Taking the Rayleigh quotient at `V^1/2 * 1` shows the slowest
mode can never decay faster than the well-mixed `eta Qs / V_total`.

Sources are piecewise constant: 06:00-22:00 two crew in Hab, one each in Lab
and Greenhouse, with the greenhouse lights on (net plant uptake 0.02 kg/h).
22:00-06:00 all four sleep in the quarters. Crew output is 1.04 kg/day
(NASA BVAD nominal), and the sleep rate is set to 0.7 of the awake rate. Each phase is
propagated exactly with the exponential of the augmented generator
`[[A dt, s dt], [0, 0]]`. One day composes to `c -> M c + b`, so the periodic
steady state is `c0 = (I - M)^-1 b`.

For a hub network the constant-source steady state has a closed form. All
CO2 leaves through the scrubber, so `c_Hab = S_total / (eta Qs)`, and
each link carries the net source behind it, so `c` rises by `flux / Q` across
it. For the sleeping crew this gives the steady-state sizing
`Q >= S_q / (c_lim - c_Node)`. The periodic sizing bisects on the peak
quarters level, found with a 1-minute scan and then golden-section refinement on the
exact map. Fan-failure and outage times use the same scan and then bisection.

ppCO2 = `rho R T / M` at 295.15 K: 1 kg/m^3 = 418.24 mmHg.

## Verification

`cargo test --release` (13 tests):

- `expm` against a rotation and a diagonal matrix. Jacobi against the
  `2 - 2cos(k pi/(n+1))` spectrum of a tridiagonal matrix.
- One-zone exponential approach and the two-zone characteristic quadratic,
  both to 1e-13.
- Hub closed form against the linear solve, for four fan sizes and both phases.
- Sum of decay rates equals `-trace(A)`. The slowest rate obeys the Rayleigh bound.
- Exact map against RK4 over each phase (1e-10). With the scrubber off, total mass rises by exactly `S t`.
- Periodic state maps to itself after one day and equals 30 days of brute-force
  RK4 from a clean cabin (1e-9). Over a day the scrubber removes exactly the
  CO2 the crew and plants add (`eta Qs <c_Hab> = <S>`).
- Crossing-time and peak searches land on the limit and on a true maximum.

`python3 verify/crosscheck.py` rebuilds the model independently in
numpy/scipy (DOP853 at rtol 1e-12, `eigvalsh`, `brentq`). It recovers the periodic state from ODE
solves rather than the exponential formula, and it compares 29 quantities
with the Rust output. All agree to 3e-12 or better, and to 4e-9 for the bisected fan size.

## Run

```
cargo run --release              # design study
cargo test --release             # verification
python3 verify/crosscheck.py     # independent scipy cross-check (needs numpy, scipy)
```

## Sample output (abridged)

```
1. Zone modes (quarters fan pair 60 m^3/h)
   time constants [h]: 7.801  0.774  0.672  0.599  0.070
   slowest 7.801 h vs. whole-habitat well-mixed V/(eta Qs) = 6.392 h

3. Periodic steady state over the 24 h schedule [mmHg]
   zone          06:00     min     max 24h mean
   Hab           1.150   1.150   1.404    1.312
   Node          1.536   1.533   1.615    1.579
   Lab           1.542   1.542   1.813    1.713
   Greenhouse    1.542   1.542   1.759    1.677
   Quarters      2.481   1.553   2.526    1.892
   mass balance: eta Qs * mean(c_Hab) = 0.160000 kg/h, mean source = 0.160000 kg/h

4. Sizing the crew-quarters fan pair for a 3.0 mmHg peak
   sleep steady-state sizing:   Q >= 37.1 m^3/h
   periodic steady-state sizing: Q >= 38.3 m^3/h  (+3 %)

5. Quarters fan pair fails at 22:00 (door gap 3 m^3/h left)
   quarters starts at 1.612 mmHg, reaches 3.0 mmHg after 55 min

6. Scrubber outage at 06:00, day schedule continuing
   Hab         reaches 3.0 mmHg after 6.45 h
   Quarters    reaches 3.0 mmHg after 7.47 h
```

## Takeaways

- The well-mixed figure (1.41 mmHg by day) hides a 2.5 mmHg sleeping
  zone. Two fan hops from the scrubber, the quarters run 1.1 mmHg above Hab.
- Steady-state sizing is slightly unconservative. The Node starts the night at
  1.61 mmHg, still carrying daytime CO2, and relaxes only on the 7.8 h
  whole-habitat mode, so the periodic answer needs 3 % more flow. The
  60 m^3/h baseline has a 0.47 mmHg margin.
- A failed quarters fan reaches the limit in under an hour, so the quarters need
  their own CO2 sensor or a fan-flow alarm. The habitat-level reading will
  not show it.
- Without the scrubber the habitat has about 6.5 h. Hab hits the limit first,
  because it holds half the daytime crew. The well-mixed estimate (6.8 h) is close
  but optimistic.
