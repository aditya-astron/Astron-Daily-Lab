# 042 - Chebyshev Ephemeris Compression for an Eccentric Orbit

**Domain:** Numerical methods / Navigation

## Problem

A ground segment has to ship a predicted ephemeris for a Molniya orbit
(a = 26 554 km, e = 0.74, 12 h period) to a user that needs position to 1 mm
and velocity to well under 1 mm/s. JPL DE files and NAIF type-2 SPK segments
do this with piecewise Chebyshev polynomials. The questions are what degree,
how to lay the segments out over a strongly eccentric orbit, and whether the
convergence rate can be predicted in advance rather than found by trial.

## Approach

**Fit and evaluation.** Each segment maps time to `x in [-1, 1]` and stores
`N+1` coefficients per axis from interpolation at first-kind Chebyshev nodes
(a type-II DCT). Evaluation is the Clenshaw recurrence. Velocity is the
derivative of the same polynomial: `d_{k-1} = d_{k+1} + 2k c_k`, scaled by
`2 / span`, so no velocity coefficients are stored.

**Truth.** An exact Kepler orbit (Newton on `E - e sin E = M`), so every error
reported is interpolation error only.

**Why the rate is predictable.** `x(t)` is entire in `E`; the only
singularities come from inverting Kepler's equation where
`dM/dE = 1 - e cos E = 0`, i.e. `E = 2 pi m +/- i acosh(1/e)`. In mean anomaly
those branch points sit at

    M* = 2 pi m +/- i eta,   eta = acosh(1/e) - sqrt(1 - e^2)  (= 0.1429 rad here)

directly above perigee. Mapping a segment onto `[-1, 1]` and taking the
Bernstein ellipse through the nearest `M*`,
`rho = |z + sqrt(z - 1) sqrt(z + 1)|`, the coefficients decay like
`k^(-3/2) rho^-k` (square-root branch point). Perigee segments see the
singularity from 0.14 rad away; apogee segments from more than pi.

**Layouts compared** at a fixed degree N, all accepted only if the max
position error on 257 points per segment is <= 1 mm:

- uniform: equal-length segments, fewest that pass;
- bisection: recursive halving of any failing segment;
- equal-rho: march from perigee, ending each segment where its predicted
  `rho` falls to `rho*`, so every segment converges at the same rate. The
  layout comes from the singularity geometry; the truth model is only used
  to pick the smallest `rho*` that passes;
- baseline: cubic Hermite through stored (r, v) at uniform knots.

## Verification (`cargo test --release`)

- Fitted coefficients of `exp(x)` and `cos(2x)` equal `2 I_k(1)` and
  `2 (-1)^k J_2k(2)` (Jacobi-Anger) to 1e-15, against A&S / scipy values.
- Clenshaw vs `cos(k acos x)`, derivative vs `k U_{k-1}(x)`, exact
  reproduction of a quintic.
- Kepler solver against Vallado Example 2-1 (E = 220.512074767522 deg),
  residuals over e up to 0.95, vis-viva along the orbit.
- Bernstein geometry: rho = 1 on the interval, minor/major axis identities,
  `eta` by hand.
- Measured envelope decay of perigee-centred segments matches the
  branch-point `rho` to within 1 %.
- The N = 13 equal-rho layout is contiguous over one period and holds 1 mm /
  1 mm/s on an independent 2003-point grid; it beats bisection, which beats
  uniform, and is more than 10x smaller than Hermite.
- Hermite error falls by ~16x per knot doubling (fourth order).

## How to run

    cargo run --release      # report
    cargo test --release     # verification

Only the Rust standard library is used.

## Sample output

```
Molniya orbit: a = 26554 km, e = 0.74, period = 717.7 min
Kepler branch points at M* = 2 pi m +/- i eta, eta = 0.1429 rad

1) Coefficient decay of x(t) on one segment, measured vs predicted rho
   centre      half-width   rho pred   rho meas
   perigee       0.15 rad      2.334      2.331  (k 8..22)
   perigee       0.30 rad      1.584      1.582  (k 14..40)
   perigee       0.60 rad      1.266      1.266  (k 26..76)
   M = 90 deg    0.60 rad      5.064      4.826  (k 6..11)
   M = 90 deg    1.00 rad      2.807      2.843  (k 6..17)
   M = 90 deg    1.50 rad      1.480      1.483  (k 16..47)
   apogee        1.00 rad      6.127      6.737  (k 6..10)
   apogee        1.50 rad      3.941      3.866  (k 6..14)
   apogee        2.50 rad      2.027      2.099  (k 9..26)

2) Doubles per orbit for 1 mm max position error (checked on 257 points per segment)
    N | uniform seg/dbl | bisection seg/dbl | equal-rho seg/dbl   rho*  | max vel err (equal-rho)
    7 |   142 /  3408    |    42 /  1008      |    30 /   720     15.55 | 0.12 mm/s
   10 |    64 /  2112    |    18 /   594      |    14 /   462      7.21 | 0.09 mm/s
   13 |    39 /  1638    |    12 /   504      |     9 /   378      4.63 | 0.09 mm/s
   18 |    22 /  1254    |    10 /   570      |     6 /   342      3.05 | 0.09 mm/s
   24 |    13 /   975    |     8 /   600      |     5 /   375      2.29 | 0.13 mm/s
   cubic Hermite on (r, v) at uniform knots: 4134 knots, 24804 doubles

3) N = 13 equal-rho layout (times from perigee)
   seg    t0 [min]   span [min]   rho    max |dr| [mm]
     0       0.00       16.24     4.63      0.213
     1      16.24       29.59     4.63      0.336
     2      45.84       67.25     4.63      0.519
     3     113.09      160.18     4.63      0.840
     4     273.27      260.03     4.63      1.000
     5     533.30      108.49     4.63      0.696
     6     641.79       46.10     4.63      0.437
     7     687.89       21.59     4.63      0.274
     8     709.47        8.24     8.27      0.000
```

Reading it: the branch-point prediction matches the measured coefficient
decay to three digits around perigee (away from perigee only a dozen
coefficients sit above round-off, so the fit is looser). At N = 13, the
DE-style degree, equal-rho segmentation needs 9 records and 378 doubles per
orbit against 39 uniform records (1638) and 4134 Hermite knots (24 804):
segments run from 16 min at perigee to 4.3 h at apogee. Above N of about 18
the extra coefficients cost more than the segments they save. Equal `rho` does
not mean equal error: the polynomial's scale grows with segment span, so the
long apogee record sets the limit. The last record is the sliver left before
the next perigee; starting the march half a segment before perigee would
absorb it.
