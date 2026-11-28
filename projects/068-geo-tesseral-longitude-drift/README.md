# 068 - Spherical-Harmonic Gravity and GEO Longitude Drift

**Domain:** Numerical methods / Orbital mechanics

## Problem

The equator is not a circle. The C22/S22 terms (plus the degree 3 and 4
tesserals) give a geostationary satellite a small along-track acceleration
that depends on its longitude. Left alone it drifts toward one of two gravity
wells and librates around it for years; kept in a slot it needs east-west
burns. Where are the wells, how long is the libration, and what does it cost
to hold a +-0.05 deg box at a given longitude?

Answering that needs a gravity model that is right to well below the 1e-8
m/s^2 signal, so half the project is the field evaluator and its checks.

C++17, header-only, no dependencies.

## Approach

**Field** (`include/gravity.hpp`), fully normalized coefficients, EGM96 to
degree and order 4:

```
U = GM/r sum_n (R/r)^n sum_m Pbar_nm(sin phi) (C_nm cos m lam + S_nm sin m lam)
```

Two acceleration paths that share nothing but the coefficients:

* `accel_spherical`: normalized Legendre functions by the sectoral seed and
  the two-term column recursion, then the spherical partials with
  `dPbar_nm/dphi = k_nm Pbar_n,m+1 - m tan(phi) Pbar_nm`,
  `k_n0 = sqrt(n(n+1)/2)`, `k_nm = sqrt((n-m)(n+m+1))`.
* `accel_cunningham`: Montenbruck & Gill V/W recursion on unnormalized
  coefficients (non-singular at the poles, overflows past degree ~150).

**Drift dynamics** (`include/geo.hpp`). A tangential acceleration changes
the semi-major axis at `da/dt = 2 a_T / n`; with `dn/da = -3n/(2a)`:

```
lam_ddot = -3 a_T(lam) / a
```

Equilibria are the zeros of `a_T`, stable where `d(lam_ddot)/d(lam) < 0`,
small-angle libration period `2 pi / sqrt(-d lam_ddot / d lam)`. `a_T` is taken
directly from `(1/r) dU/dlam`; projecting the Cartesian vector loses about
five digits to the central term.

The 1D model is checked against a full 3D orbit: RK4 in inertial space with
the Earth rotating uniformly, started on a circular orbit whose radius
balances the local radial gravity (all harmonics). Longitude is the
sidereal-day running mean, so the daily eccentricity wobble drops out.

**Station keeping.** A cycle is a parabola in longitude: leaving the downhill
edge at `v = sqrt(2 A W)` against drift acceleration `A` carries the
satellite exactly to the uphill edge and back in `2v/A`. The reversal costs
`dv = 2 v a / 3`, so

```
dv rate = a A / 3 = |a_T|          (independent of the box width)
```

The simulator uses paired tangential burns 12 h apart (no net eccentricity to
first order) and closes the loop on the drift estimated from the daily mean:
a fixed open-loop reversal lets each cycle's timing error feed the next and
the box walks out to several tenths of a degree.

Scope: Earth gravity only. Sun/Moon (inclination drift, north-south keeping)
and solar pressure (eccentricity) are left out.

## Verification (`tests/verify.cpp`)

| Check | Result |
|---|---|
| J2, J3, J4, lambda22 from EGM96 vs 1.08263e-3, -2.5327e-6, -1.6196e-6, -14.93 deg | match |
| Pbar_20..Pbar_40 vs closed forms | 9e-16 |
| Orthonormality to degree 180, 200-point Gauss-Legendre | 5e-14 |
| EGM96 4x4: spherical vs Cunningham / vs 4th-order FD of U | 8e-16 / 4e-12 |
| Random degree-60 field: spherical vs Cunningham / vs FD of U | 5e-15 / 3e-11 |
| J2-only field vs textbook closed form | 6e-16 |
| Degree 2: equilibria at lam22 + k 90 deg, curvature `6 sqrt15 GM R^2 J22 / r^5` | 3e-9 deg, 2e-10 |
| EGM96 4x4 equilibria vs commonly quoted 75.1 E, 161.9 E, 105.3 W, 11.5 W | within 0.21 deg |
| Small-angle period vs 1D integration; 1D vs 3D orbit at 20 deg | < 1e-4, 2e-5 |
| 3D station keeping at 30 E: dv rate vs `|a_T|` x 1 yr; cycle vs `2 sqrt(2W/A)`; box | 0.01 %, 2.3 %, 0.048 deg |

The 1D and 3D periods differ by a steady 0.002 %, the size of the J2 term
`J2 (R/a)^2 ~ 2.5e-5` that the Keplerian `dn/da` leaves out.

## Results

* Stable wells at **74.94 E** and **105.09 W**, unstable points at 161.90 E and
  11.52 W. Degree 2 alone puts them exactly 90 deg apart on the axes of the
  equatorial ellipse; J3x moves the unstable ones by 3 deg, degree 4 by
  another 0.1-0.2 deg.
* The wells are not equally deep: curvature -7.24e-5 vs -4.80e-5 per day^2,
  so small librations take **738 days** about 75 E and **907 days** about
  105 W. A 60 deg release stretches both to ~1145 days and the swing becomes
  lopsided (-65.8/+60 deg about 75 E).
* Peak drift acceleration 0.002 deg/day^2 near 117 E. East-west keeping costs
  `|a_T|` per unit time: 0.67 m/s/yr at 0 E, 1.84 at 30 E, **2.07 m/s/yr** at
  the worst slot, and nothing at the four equilibria. The 3D box-keeping runs
  reproduce the rate to 0.1 % and hold the daily mean inside +-0.05 deg.

## Run

```
make            # builds geo_longitude and verify
./verify        # ~6 s, exits non-zero on any failure
./geo_longitude # ~35 s, full study; output in sample_output.txt
```

## Sample output (excerpt)

```
Equilibria, field to degree 4:
    74.939 E  ( 74.939 E)  stable    d(lam_ddot)/d(lam) = -7.239e-05 1/day^2
   161.904 E  (161.904 E)  unstable  d(lam_ddot)/d(lam) = +6.515e-05 1/day^2
   254.906 E  (105.094 W)  stable    d(lam_ddot)/d(lam) = -4.796e-05 1/day^2
   348.482 E  ( 11.518 W)  unstable  d(lam_ddot)/d(lam) = +5.522e-05 1/day^2

Free libration, released at rest (period in days):
  well      release      1D model    3D orbit    diff
    74.94 E  small-angle      738.5
            +  20 deg        767.7       767.7   -0.002 %   (swing -20.3/+20.0 deg)

East-west box keeping, +-0.05 deg, paired tangential burns, 3 years:
  slot    predicted |a_T|*yr   simulated     cycle     burns   max excursion
     30 E      1.843 m/s/yr       1.842      20.7 d     106     0.0480 deg
     60 E      1.041 m/s/yr       1.041      27.7 d      80     0.0489 deg
```

## References

* O. Montenbruck, E. Gill, *Satellite Orbits*, Springer 2000, sec. 3.2.
* F. Lemoine et al., *The Development of the Joint NASA GSFC and NIMA
  Geopotential Model EGM96*, NASA/TP-1998-206861.
* E. M. Soop, *Handbook of Geostationary Orbits*, Kluwer 1994.
