# 009 - Multiplicative EKF Fusing a Drifting Gyro with a Star Tracker

**Domain:** Kalman filtering and sensor fusion

## Problem

A spacecraft slewing at about 0.14 deg/s carries a navigation-grade gyro
(ARW 3.2e-7 rad/s^0.5, RRW 3.2e-10 rad/s^1.5, initial bias up to 0.15 deg/hr)
sampled at 10 Hz and a 6 arcsec star tracker at 1 Hz. Estimate attitude and
gyro bias, show the filter is statistically consistent, and check its
steady-state accuracy against the analytic solution.

## Approach

The MEKF keeps a global unit quaternion and a 6-state error vector
`x = [dtheta, dbeta]` (body-frame small rotation, bias error), so the covariance
is 6x6 and never singular.

Gyro model (Farrenkopf):

```
w_meas = w + beta + eta_v,     beta_dot = eta_u
```

Error dynamics with `w_hat = w_meas - beta_hat`:

```
dtheta_dot = -[w_hat x] dtheta - dbeta - eta_v
dbeta_dot  = eta_u
```

- Propagation: `q <- q (x) exp(w_hat dt / 2)`; `Phi` and `Qd` from Van Loan's
  matrix exponential, which at zero rate reduces to the closed form
  `Q11 = (sv^2 dt + su^2 dt^3/3) I`, `Q12 = -su^2 dt^2/2 I`, `Q22 = su^2 dt I`.
- Update: innovation `y = 2 vec(q_est* (x) q_meas)`, `H = [I 0]`, Joseph-form
  covariance update, then reset by folding `dtheta` into the quaternion.
- Truth gyro is rate-integrating: the mean bias over each step plus noise of
  variance `sv^2/dt + su^2 dt/12`, so the integrated angle has exactly the
  statistics the filter assumes.
- Analytic reference: per axis at zero rate, the 10 gyro steps between star
  tracker updates combine into one exact 1 s step; the discrete algebraic
  Riccati equation then gives the steady-state attitude and bias variance.

## Verification (`test_mekf.py`)

- Quaternion error round-trip to 1e-12.
- Van Loan `Phi`, `Qd` equal the Farrenkopf closed form at zero rate; at
  non-zero rate the attitude block equals the Rodrigues rotation transpose.
- DARE steady state equals a brute-force Riccati iteration with 10 substeps.
- In the `su -> 0` limit the result equals the scalar random-walk closed form
  `P- = q/2 + sqrt(q^2/4 + qR)`, `P+ = P- R/(P- + R)`.
- 12-run Monte Carlo: ensemble NEES inside the 99% chi-square band, final
  sigmas within 2% of analytic, normalised attitude and bias RMS near 1,
  every run's bias learnt to within 3 sigma.

## Run

```
pip install numpy scipy pytest
python3 mekf.py
python3 -m pytest -q
```

## Sample output

```
MEKF: 10 Hz gyro + 1 Hz star tracker, 1 h, slewing at ~0.14 deg/s
  gyro ARW 3.162e-07 rad/s^0.5, RRW 3.162e-10 rad/s^1.5, ST 6.0 arcsec
  steady state (t > 1200 s)       x        y        z    analytic
  attitude RMS [arcsec]        0.777    0.549    0.603     0.651
  filter 1-sigma [arcsec]      0.650    0.651    0.651
  bias RMS [deg/hr]           0.0017   0.0017   0.0023    0.0022
  filter 1-sigma [deg/hr]     0.0022   0.0022   0.0022
  mean NEES (6 dof)         5.35
  attitude improvement over raw star tracker: 9.22x
```

The gyro carries the attitude between star fixes so the fused estimate is
about 9x better than a single star tracker measurement, and the bias is
pinned to ~0.002 deg/hr, roughly 50x below its initial value.
