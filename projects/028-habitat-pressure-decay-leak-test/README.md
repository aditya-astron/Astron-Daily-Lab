# 028 - Habitat Pressure-Decay Leak Test

**Domain:** Analog habitat systems

## Problem

A sealed analog habitat module (300 m^3, rigid hull) is pumped to +3 kPa,
sealed, and logged for 72 h. Find its leak rate at the 1 kPa operating
differential, the leak's flow exponent, and the equivalent hole size, with an
honest uncertainty.

The catch is the weather. At 101 kPa and 293 K one kelvin of air temperature
moves the hull pressure by p/T = 346 Pa, and the barometer moves by hundreds
of Pa a day. In the reference log a perfectly sealed hull would still show a
1100 Pa swing in differential pressure, while the leak itself takes about
1800 Pa off over three days. Reading the leak off the pressure trace alone is
biased by 17 % here, even when the flow exponent is given.

## Approach

Track the gas inventory instead of the pressure (ideal gas, rigid hull):

```
m = (p_e + dp) V / (R T)
```

and model the leak with the envelope power law (ASTM E779 form),

```
Q  = C |dp|^n sign(dp)          n = 0.5 orifice ... 1.0 laminar crack
dm/dt = -rho_up Q               rho_up = m/V outflow, p_e/(R T_e) inflow
```

The fit is output-error: integrate that ODE with RK4 through the logged
temperature and barometer, and adjust (ln C, n, m0) by Levenberg-Marquardt
until the model inventory matches the measured one. Parameter covariance is
s^2 (J^T J)^-1; the uncertainty on Q(dp_ref) follows by the delta method since
ln Q = ln C + n ln dp_ref. The equivalent orifice is A = Q / (Cd sqrt(2 dp/rho))
with Cd = 0.61.

Two baselines run on the same log: the uncorrected pressure-trace fit
(isothermal decay, `dp^(1-n)` linear in t) and the textbook two-point
temperature-compensated method over whole days, which gets the mass loss right
but gives a single rate with no exponent to extrapolate with.

Expected residual, by hand: temperature noise 0.02 K gives
371.7 kg x 0.02/293 = 25.4 g, barometer 2 Pa gives 7.3 g, transducer 0.5 Pa
gives 1.8 g; RSS 26.5 g. The fit reports 26.3 g.

## Verification (`leak_test.go`)

1. Orifice area and flow law against hand calculations.
2. Sealed hull (C = 0): dp swings 1116 Pa with the weather, inventory stays
   constant to 1e-13 relative. This is the premise of the method.
3. Integrator vs. an independent Simpson quadrature of the decay time,
   `t = int V / (C (p_e + x) x^n) dx`, for n = 0.62 (error 1e-13), and vs. the
   closed form for n = 1.
4. RK4 order: error ratio 16.5 on step halving under time-varying weather.
5. Inflow branch uses exterior density.
6. Noiseless log: n and Q(1 kPa) recovered to a few parts per million.
7. Same noisy log: pressure-only error > 5 %, inventory fit < 1 %.
8. Monte Carlo, 150 noise realisations: bias of Q and n within 4 standard
   errors, empirical scatter within 25 % of the reported sigma
   (Q: 0.244 % actual vs 0.254 % predicted; n: 0.0037 vs 0.0038),
   n inside 2 sigma in 144/150 runs.

## How to run

Requires Go 1.24+, no external modules.

```
go test -v ./...          # add -short to skip the Monte Carlo
go run ./cmd/leaktest
```

## Sample output

```
Habitat pressure-decay leak test (reference scenario)
  volume 300 m^3, sealed at +3000 Pa, 72 h log every 300 s
  dp: 3000 Pa -> 1194 Pa, inventory 371.72 kg -> 363.86 kg
  temperature sensitivity p/T = 346 Pa/K; barometer front -400 Pa

Truth: C = 2.300e-07 m^3/(s Pa^n), n = 0.620
  Q(1000 Pa) = 1.666e-05 m^3/s = 1.733 kg/day = 0.466 %/day of inventory

1) Pressure trace only (no T or barometer correction, n given):
   Q(1000 Pa) = 1.383e-05 m^3/s  error -17.0 %

2) Two-point temperature-compensated method, 24 h windows:
   day 1: mdot = 3.1067 kg/day at mean dp 2478 Pa  (true mass loss 3.1043, error +0.1 %)
   day 2: mdot = 2.6179 kg/day at mean dp 1888 Pa  (true mass loss 2.6073, error +0.4 %)
   day 3: mdot = 2.0623 kg/day at mean dp 1411 Pa  (true mass loss 2.1460, error -3.9 %)
   (one rate per window: no exponent, so no extrapolation to other dp)

3) Inventory fit (Levenberg-Marquardt, 7 iterations, residual 26.3 g):
   n = 0.621 +/- 0.004   (truth 0.620)
   C = 2.282e-07  (truth 2.300e-07)
   Q(1000 Pa) = 1.6658e-05 +/- 4.2e-08 m^3/s  error -0.02 % (-0.1 sigma)
   leak = 1.733 kg/day = 0.466 %/day at 1000 Pa
   equivalent orifice (Cd 0.61) at 1000 Pa: 0.670 mm^2, d = 0.924 mm
   equivalent orifice (Cd 0.61) at 3000 Pa: 0.765 mm^2, d = 0.987 mm
```

## Limits

The logged temperature is assumed to be the volume-averaged gas temperature.
A stratified module with sensors in the wrong places produces a slow,
correlated temperature error that the fit cannot separate from leakage; the
usual remedies are more sensors, a fan to mix the air, and windows that span
whole days. Hull compliance (volume change with pressure) is neglected, which
is fine for a rigid module but not for an inflatable one.
