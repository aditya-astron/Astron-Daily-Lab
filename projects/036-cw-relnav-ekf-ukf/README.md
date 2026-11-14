# 036 - Rendezvous Relative Navigation: EKF vs UKF vs the Cramer-Rao Bound

**Domain:** Kalman filtering and sensor fusion

## Problem

A chaser hands over from relative GPS to its own sensors 8 km behind and
500 m below a target in a 420 km circular orbit, on a coelliptic drift that
brings it to about 1 km in 1.5 orbits. The hand-over state is uncertain to
100 m radial, 500 m along-track and 5 cm/s. Two sensor suites:

- **Lidar**: range (1 m) plus azimuth and elevation (1 mrad), every 10 s.
- **Camera only**: the two angles, nothing else.

Questions: are an EKF and a UKF statistically consistent and as accurate as
any estimator can be? And with angles only, can the chaser ever learn its
range?

## Approach

State `[x y z vx vy vz]` in the target's Hill frame (x radial, y along-track,
z orbit normal), Clohessy-Wiltshire dynamics with white acceleration noise
of PSD `q = 1e-10 m^2/s^3` (differential drag, small leaks):

    x'' = 3n^2 x + 2n y' + w,   y'' = -2n x' + w,   z'' = -n^2 z + w

- Exact discrete model: closed-form CW transition matrix `Phi(dt)` and
  Van Loan's `Q = int Phi G q G^T Phi^T`.
- Line of sight `u = -r`; `az = atan2(u_x, u_y)`, `el = asin(u_z/|u|)`,
  analytic Jacobian, angle residuals wrapped.
- **EKF** with Joseph-form update. **UKF** with the scaled unscented
  transform (alpha 1, beta 2, kappa 0); angle means are taken on residuals
  about the central sigma point. Dynamics are linear, so both share the
  exact prediction and differ only in the measurement update.
- **RTS smoother** for post-pass reconstruction of the hand-over state.
- **Cramer-Rao bound** along each true trajectory, i.e. the information
  filter run on the truth:
  `J+ = (Q + Phi J^-1 Phi^T)^-1 + H(x_true)^T R^-1 H(x_true)`.
  The Bayesian PCRB (Tichavsky 1998), which averages `H^T R^-1 H` over the
  ensemble first, is also printed; it is looser by Jensen's inequality.
- **Angles-only observability.** Angles are homogeneous of degree zero in
  `r` and CW is linear, so `alpha * x(t)` yields the same angles for any
  `alpha > 0`. Linearised, `H(x_k) Phi_k x0 = H(x_k) x_k = 0`: the initial
  state itself is the null space of the observability Gramian, and the only
  range information on a coast is the prior. A known impulse adds a term
  that does not scale, which removes the null space.

## Checks (`test_relnav.py`)

- CW `Phi` equals `expm(A t)` (5e-14), is a semigroup, and matches direct
  ODE integration to 6e-10 m.
- Van Loan `Q` equals adaptive quadrature (2e-16) and the double-integrator
  limit `q dt^3/3, q dt^2/2, q dt` to `O((n dt)^2)`.
- Lidar Jacobian against central differences (1e-9).
- Angles-only coast: scaled trajectory gives identical angles (1e-16 rad);
  `|W x0| / (lambda_max |x0|) = 2e-16`, next eigenvalue 1e-4. A 0.1 m/s
  radial burn makes the Gramian full rank.
- UKF reproduces the linear Kalman filter exactly on a linear position
  sensor (8e-12); the Riccati recursion converges to `solve_discrete_are`.
- With `Q = 0`, the RTS-smoothed initial state and covariance equal the
  batch MAP solution (6e-11).
- Monte Carlo, lidar: run-averaged NEES and NIS inside 99 % chi-square
  bounds on more than 90 % of steps; MSE within 1 % of the CRB.
- Angles only: with a 5x looser prior the coast range bound grows 4.6x
  (208 to 964 m); after the burn it barely moves (81 to 88 m).

## Run

    pip install numpy scipy
    python3 relnav.py        # ~1 min, 100-run Monte Carlo for each case
    python3 test_relnav.py   # ~12 s

## Sample output

```
Rendezvous relative navigation, target at 420 km circular
  n = 1.126378e-03 rad/s, T = 5578.2 s, dt = 10 s, 837 steps (1.50 orbits)
  nominal coelliptic drift: r0 = [ -500. -8000.   100.] m -> r_end = [-500.  -929.2 -100. ] m

Case A: lidar (1 m, 1 mrad), q = 1e-10 m^2/s^3, 100 Monte Carlo runs
  99% bounds on run-averaged NEES [5.15, 6.93], NIS [2.41, 3.67]
  filter  NEES in  NIS in  mean NEES  mean NIS  pos RMSE end  vel RMSE end  MSE/bound
  EKF      99.4%   98.2%       6.07      3.00        0.52 m     0.89 mm/s     1.063
  UKF     100.0%   98.3%       6.00      2.99        0.52 m     0.89 mm/s     1.063
  per-trajectory CRB at end: pos 0.51 m, vel 0.88 mm/s (MSE/bound averaged over the last 2/3 of the arc)
  Bayesian PCRB at end: pos 0.40 m (looser: averages information over the range spread)
  RTS smoother, one run: hand-over state pos sigma 11.3 m filtered -> 0.61 m smoothed (actual error 0.95 m)

Case B: camera angles only (1 mrad), 100 Monte Carlo runs
  coast:
    t = T/2 + 100 s  CRB range sigma  215.4 m (5x prior: 1066.9 m) | EKF RMSE  267.3 m | UKF RMSE  257.7 m
    end              CRB range sigma  189.0 m (5x prior:  894.5 m) | EKF RMSE  334.6 m | UKF RMSE  322.6 m
    mean NEES (ideal 6): EKF 6.87, UKF 6.52
  radial burn 0.10 m/s at t = 2789 s:
    t = T/2 + 100 s  CRB range sigma  120.5 m (5x prior:  142.9 m) | EKF RMSE  154.4 m | UKF RMSE  145.0 m
    end              CRB range sigma   72.1 m (5x prior:   78.2 m) | EKF RMSE  156.6 m | UKF RMSE  139.8 m
    mean NEES (ideal 6): EKF 6.50, UKF 6.25
```

With lidar the problem is nearly linear: once the first few range returns
arrive, both filters sit on the bound and the UKF buys nothing. With angles
only, the coast range error is set entirely by the hand-over prior (loosen
it 5x and the bound loosens 5x), and both filters run 1.2-2.2x above the
bound and slightly overconfident because they linearise about the wrong
range. A 10 cm/s radial burn makes range data-limited, and the UKF handles
the larger post-burn nonlinearity better (140 m vs 157 m).

**Limits.** Circular-orbit CW model (no J2, eccentricity or differential
drag model beyond white noise), known attitude, no sensor biases or
outliers, burn executed perfectly. Rectilinear CW coordinates misplace the
chaser radially by about `y^2 / 2a`, roughly 5 m at 8 km, comparable to one
angle sample there; curvilinear coordinates or a full nonlinear model are
needed for longer arcs.
