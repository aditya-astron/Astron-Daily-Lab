# 026 - Batch Least-Squares Orbit Determination from Range and Range-Rate

**Domain:** Navigation

## Problem

A 500 km, 51.6 deg LEO spacecraft is tracked for 24 h by three ground
stations (Bengaluru, Hartebeesthoek, Wallops) above a 10 deg elevation mask.
Each pass gives two-way range (sigma 10 m) and range-rate (sigma 1 mm/s) every
10 s. Starting from an initial guess that is ~2.7 km and ~2.7 m/s off,
estimate the epoch state and a covariance that can actually be trusted.

## Approach

Dynamics are two-body + J2 in ECI. The state transition matrix is integrated
with the state (DOP853, rtol 1e-11):

    dPhi/dt = A Phi,   A = [[0, I], [G, 0]],   G = da/dr  (analytic, incl. J2)

Stations are placed on the WGS-84 ellipsoid and spun at the Earth rate.
With `d = r - R`, `w = v - Vs`:

    rho    = |d|                  d rho/dr    = d/rho
    rhodot = d.w / rho            d rhodot/dr = w/rho - (d.w) d/rho^3
                                  d rhodot/dv = d/rho

Each measurement row is mapped to the epoch with `H = H~ Phi(t, t0)`, and the
batch solves the Gauss-Newton normal equations

    (H^T W H) dx0 = H^T W dy,    P0 = (H^T W H)^-1

with the velocity columns scaled by 1e3 so the normal matrix stays well
conditioned. It iterates until the update is below 1 mm and 1 um/s.

**Limits.** No light time, troposphere, clock bias or precession/nutation, and
the truth comes from the same force model, so the result shows the
estimator's noise-limited accuracy, not model error.

## Verification

`test_orbit_det.py` checks:

- WGS-84 equator and pole points (polar radius 6356.752 km)
- equatorial J2 gravity against `mu/r^2 (1 + 1.5 J2 (Re/r)^2)`
- analytic gravity gradient against central differences, and its symmetry
- range/range-rate partials against central differences
- the STM after 3 h against perturbed propagations, and `det(Phi) = 1`
- noise-free data: the truth is recovered to < 1 mm and < 1 um/s from a
  2.7 km / 2.7 m/s guess
- 60-trial Monte Carlo: summed NEES inside the 99% chi-square(360) interval,
  per-axis sample variance against the formal variance, mean post-fit
  reduced chi-square within 0.02 of 1

## Run

    pip install numpy scipy
    python3 orbit_det.py
    python3 test_orbit_det.py      # or: pytest -q  (about 1 min)

## Sample output

```
Arc 24 h, 10 s sampling, 448 range/range-rate pairs
  Bengaluru         91 pairs
  Hartebeesthoek   148 pairs
  Wallops          209 pairs

Gauss-Newton iterations:
  iter 1: cost 143462206928808.2188  |dr|    8311.1673 m  |dv|   13590.1396 mm/s
  iter 2: cost 235831359442.0198  |dr|    9862.9172 m  |dv|   11161.4136 mm/s
  iter 3: cost 71769438674.8831  |dr|      20.6445 m  |dv|      21.7142 mm/s
  iter 4: cost     48001.2810  |dr|       5.1348 m  |dv|       5.8023 mm/s
  iter 5: cost       920.4465  |dr|       0.0000 m  |dv|       0.0000 mm/s

Post-fit residual RMS: range 10.23 m (sigma 10 m), range-rate 1.004 mm/s (sigma 1 mm/s)
Weighted cost 920.4 for 890 dof (reduced chi-square 1.034)

Epoch state error vs formal 1-sigma (RTN):
  R  pos   -0.012 m  (1s  0.009 m)   vel   0.0285 mm/s  (1s 0.0472 mm/s)
  T  pos   -0.031 m  (1s  0.044 m)   vel   0.0134 mm/s  (1s 0.0099 mm/s)
  N  pos    0.028 m  (1s  0.039 m)   vel  -0.0184 mm/s  (1s 0.0416 mm/s)
  NEES (6 dof, expect ~6): 5.41
```

The first step overshoots (the 24 h arc is far from linear at a few km of
along-track error), yet Gauss-Newton settles in five iterations. Range-rate
at 1 mm/s dominates the information: the epoch position comes out at the
centimetre level despite 10 m range noise.
