# 006 - Optimal Mass Staging of a Three-Stage Launcher

**Domain:** Rocket propulsion and staging

## Problem

Put a 10 t payload into a 200 km LEO with a three-stage vehicle: a LOX/RP-1
booster, a LOX/RP-1 second stage and a LOX/LH2 upper stage. The design
velocity budget is 9.0 km/s (about 7.79 km/s circular speed plus 1.6 km/s of
gravity, drag and steering losses, less the eastward launch-site rotation
credit). Given each stage's Isp and structural ratio, how should the 9 km/s be
split between stages to minimise lift-off mass, and how much does that buy
over a naive equal split?

## Approach

Restricted staging model (Curtis, *Orbital Mechanics for Engineering
Students*, sec. 11.6). For stage `i`, with everything above it as its payload:

```
eps_i = m_inert / (m_inert + m_prop)        structural ratio
n_i   = m_ignition / m_burnout              mass ratio
dv    = sum c_i ln n_i,  c_i = g0 Isp_i
m_ign_i / m_above_i = n_i (1 - eps_i) / (1 - eps_i n_i)
```

Minimising `ln(m0/m_PL)` subject to the dv constraint with a Lagrange
multiplier `eta` gives the optimality condition

```
1 / (1 - eps_i n_i) = eta c_i   =>   n_i = (eta c_i - 1) / (eta c_i eps_i)
```

`eta` is the root of `sum c_i ln n_i(eta) = dv`. The left side is strictly
increasing in `eta`, so the solver brackets it between
`max 1/(c_i (1 - eps_i))` (some `n_i -> 1`) and infinity (`n_i -> 1/eps_i`, the
zero-payload limit) and bisects. Stage masses are then built top-down from
the payload. Inputs that need more than the zero-payload dv, or whose optimum
would drive a stage to `n <= 1` (the stage should be deleted), are reported
as errors instead of returning nonsense.

The model sizes mass only. It does not see thrust-to-weight or how gravity
loss depends on stage-1 burn time, which is why it puts the least dv on the
low-Isp booster. Real boosters carry more dv than this for thrust reasons.

## Verification

`test_staging` checks:

- Single stage by hand: `c = 3000 m/s`, `eps = 0.1`, `dv = c ln 4` gives
  `n = 4`, `m0 = 6 m_PL`, 4500 kg propellant and 500 kg inert per tonne.
- 2, 3 and 4 identical stages against the closed form
  `n = exp(dv/(N c))`, `m0/m_PL = [n(1-eps)/(1-eps n)]^N`.
- Stage dv sums to the target; inert + propellant + upper mass rebuilds every
  ignition mass; the rocket equation recomputed from masses returns each
  stage dv.
- The first-order condition holds per stage.
- A 5 m/s grid search over the dv split never beats the optimum and lands
  within 0.01% of it.
- Out-of-range dv and mass ratios are rejected.

The case-study optimum was also cross-checked with an independent
Nelder-Mead minimisation in SciPy (same split to 0.1 m/s, same 182220 kg).

## Build and run

```
make test    # build and run the checks
make run     # case study
```

Needs a C++17 compiler, no other dependencies.

## Sample output

```
Optimal staging: payload 10000 kg, dv 9000 m/s (stack limit 26487 m/s)

Lagrange optimum
  stage         Isp   eps       n  dv [m/s]    ign [kg]   prop [kg]  inert [kg]
  S1 RP-1       300  0.05  1.4631    1119.5      182220       57673        3035
  S2 RP-1       345  0.07  2.7721    3449.6      121511       77677        5847
  S3 LH2        450  0.14  2.7293    4430.8       37987       24069        3918
  total dv 9000.0 m/s | lift-off 182220 kg | payload fraction 0.0549
  eta = 3.667332e-04 s/m

Equal dv split (baseline)
  stage         Isp   eps       n  dv [m/s]    ign [kg]   prop [kg]  inert [kg]
  S1 RP-1       300  0.05  2.7724    3000.0      194986      124655        6561
  S2 RP-1       345  0.07  2.4271    3000.0       63770       37496        2822
  S3 LH2        450  0.14  1.9735    3000.0       23452       11568        1883
  total dv 9000.0 m/s | lift-off 194986 kg | payload fraction 0.0513

  optimum saves 12767 kg of lift-off mass (6.5%)

Sensitivity, eps + 0.01 per stage (re-optimised):
  S1 RP-1    lift-off    183086 kg (+0.5%)
  S2 RP-1    lift-off    186272 kg (+2.2%)
  S3 LH2     lift-off    188007 kg (+3.2%)

Stage count, identical stages (Isp 350 s, eps 0.08):
  N=1  infeasible (required dv exceeds the zero-payload limit of the stack)
  N=2  lift-off       235616 kg  payload fraction 0.0424
  N=3  lift-off       202985 kg  payload fraction 0.0493
  N=4  lift-off       192591 kg  payload fraction 0.0519
  N=5  lift-off       187514 kg  payload fraction 0.0533
  N=6  lift-off       184511 kg  payload fraction 0.0542
```
