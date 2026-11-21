# 053 - Li-ion state-of-charge EKF for a LEO battery: sensor bias, correlated noise and the LFP plateau

**Domain:** Kalman filtering and sensor fusion / Spacecraft power systems

## Problem

A LEO battery cycles every 95 minutes, and the power subsystem needs state of
charge to decide load shedding and end-of-charge. Coulomb counting integrates
the current sensor, so a 20 mA Hall-sensor offset on a 3 Ah cell walks SOC by
1 % per orbit and never recovers. The terminal voltage pins SOC through the
OCV curve, but only as well as the curve's slope allows. This project builds
an EKF that fuses both, estimates the sensor bias, handles the noise
correlation the shared current sensor creates, and compares a sloped NMC cell
with a flat-plateau LFP cell.

## Approach

**Cell model.** Second-order Thevenin, I > 0 discharging:

```
dz/dt   = -I / Q
dv_i/dt = -v_i / (R_i C_i) + I / C_i        (i = 1, 2; tau = 30 s, 1000 s)
V       = OCV(z) - v1 - v2 - R0 I
```

With I held over a sample (ZOH) the RC branches discretise exactly:
`v+ = a v + R (1 - a) I`, `a = exp(-dt / RC)`.

**Sensors.** `Im = I + b + n` (bias b, white n of 10 mA), `y = V + e` (2 mV).
The filter state is `x = [z, v1, v2, b]`, driven by `Im - b`.

**Correlated noise.** The same n enters the dynamics (`w = Bn n`,
`Bn = [dt/Q, -R1(1-a1), -R2(1-a2), 0]`) and the voltage residual
(`v = R0 n + e`), so `S = E[w v] = Bn R0 sigma_n^2`. Writing
`w = S v / R + w~` gives an uncorrelated system with

```
x+ = f(x, Im) + S (y - h(x)) / R
F~ = F - S H / R,   Q~ = Q - S S^T / R
```

applied after the usual measurement update, with a Joseph-form covariance.

**OCV curves.** NMC: polynomial plus log knees, 11 to 21 mV per % SOC.
LFP: 0.5 mV per % on the plateau with exponential knees, so one 2 mV sample
resolves SOC to only about 4 %.

**Scenario.** 15 orbits, 35 min eclipse at 1.4 A with a +/-0.4 A payload
square wave, 1.0 A charge tapering into 90.5 %. 100 Monte Carlo runs per
chemistry with a 10 % (1 sigma) initial SOC guess error and a N(0, 20 mA)
sensor bias. Truth and filter share the cell parameters, so the errors
below are a floor set by the sensors and the estimator, not by model error.

## Verification (`make test`)

1. ZOH update vs. RK4 with 20000 substeps, dt from 0.1 to 60 s: 1e-16 V.
2. Analytic dOCV/dz vs. central differences to 1e-8; OCV monotone; bisection
   inverse round trip to 1e-14.
3. Linear OCV: the filter's x_N and P_N match the brute-force Gaussian
   conditional of x_N on y_0..y_{N-1} (40 steps, joint covariance over 124
   noise terms, Cholesky solve) to 5e-13 sigma and 2e-12. Dropping the cross
   term misses the exact covariance by about 9 (correlation units), which is
   what makes the check meaningful.
4. Coulomb-count drift equals `-b t / Q` to 1e-12 over 10 orbits.
5. Linear-OCV Monte Carlo (400 runs): mean NEES at five instants inside the
   99 % chi-square band for n = 4.
6. Short LEO Monte Carlo reproduces the qualitative findings below.

Hand check on the coulomb counter: 15 orbits x 5700 s x 20 mA / 10800 As =
15.8 % drift, combined with the 10 % start error gives sqrt(15.8^2 + 10^2) =
18.7 % expected RMS; the Monte Carlo gives 18.0 %.

## Results

- **NMC.** The bias-estimating EKF holds SOC to 0.003 % RMS with NEES 3.6
  (n = 4) and estimates the 20 mA offset to 0.07 mA. Coulomb counting is off
  by 18 %.
- **No bias state is worse than it looks.** Inflating the current noise to
  cover the 20 mA spec does not help: a constant bias is not white noise, the
  filter is 16x (NMC) to 50x (LFP) overconfident and only 14 % / 5 % of samples fall
  inside its 3-sigma bound.
- **LFP and the knee trap.** Starting from a guess on the plateau while the
  cell sits at the 90 % knee, 21 of 100 runs push the SOC innovation into the
  bias (H_b = R0 = 30 mV/A dwarfs the plateau slope), briefly estimating
  -2.9 A of offset. They recover SOC but keep a bias error of ~1.5 mA against
  a 0.12 mA sigma, so the filter stays overconfident (mean NEES 56). Clamping
  the bias to the sensor spec was tried and made it worse: the residual then
  drives SOC far out of range.
- **Fix: initialise from the first sample.** With the RC branches relaxed,
  `z0 = OCV^-1(y0 + R0 Im0)` with a sigma from the noise and local slope
  (capped at the prior). LFP then gives 0.009 % RMS, NEES 3.5, no outliers.
  This works because the pack starts at the knee; started mid-plateau the
  inversion has nothing to offer and the sigma cap applies, so the guess-init
  risk returns.
- The correlated-noise term is exact but small here (R0 sigma_n = 0.3 mV vs.
  2 mV ADC noise); in the LEO runs it moved LFP RMS from 0.051 % to 0.044 %.
  It matters with a noisy current sensor (check 3 uses 0.2 A).

Not modelled: OCV hysteresis (significant on LFP), temperature, capacity
fade and parameter mismatch.

## How to run

    make run     # report, about 18 s
    make test    # verification, about 4 s

C++17, standard library only.

## Sample output

```
Cell: Q = 3.0 Ah, R0 = 30 mOhm, tau1 = 30 s, tau2 = 1000 s
Sensors: sigma_V = 2.0 mV, sigma_I = 10 mA white, bias ~ N(0, 20 mA)

1) Exact ZOH vs RK4 (100 substeps), 3 LEO orbits: max |dv| = 2.59e-15 V

2) OCV slope and the SOC resolution of one voltage sample (sigma_V / slope)
     SOC  |   NMC mV/%   1-sample %  |   LFP mV/%   1-sample %
      10% |      14.80       0.135   |      11.52       0.174
      30% |      10.79       0.185   |       0.89       2.240
      50% |      11.15       0.179   |       0.51       3.889
      70% |      14.20       0.141   |       0.53       3.787
      90% |      21.11       0.095   |       4.60       0.434

3) LEO cycling, 15 orbits (23.75 h), 100 runs per chemistry
   initial SOC guess error 10 % (1 sigma); statistics over the last orbit

   NMC
   method                 RMS err %   max err %   mean sigma %   in 3 sigma
   coulomb count             17.987      48.939         -            -
   EKF, no bias state         0.063       0.212        0.004        14.3 %
   EKF + bias                 0.003       0.012        0.003        99.7 %
   EKF + bias, OCV init       0.003       0.012        0.003        99.7 %
   EKF + bias           : NEES mean 3.64 / per-run median 3.23 (n = 4), runs > 3n: 0,
                          worst run RMS 0.006 %, bias RMS 0.07 mA, median time to
                          stay within 1 % SOC: 0 s
   EKF + bias, OCV init : NEES mean 3.64 / per-run median 3.23 (n = 4), runs > 3n: 0,
                          worst run RMS 0.006 %, bias RMS 0.07 mA, median time to
                          stay within 1 % SOC: 0 s

   LFP
   method                 RMS err %   max err %   mean sigma %   in 3 sigma
   coulomb count             17.987      48.939         -            -
   EKF, no bias state         0.538       1.794        0.011         4.9 %
   EKF + bias                 0.044       0.190        0.010        82.9 %
   EKF + bias, OCV init       0.009       0.030        0.010       100.0 %
   EKF + bias           : NEES mean 55.64 / per-run median 3.56 (n = 4), runs > 3n: 21,
                          worst run RMS 0.154 %, bias RMS 0.62 mA, median time to
                          stay within 1 % SOC: 500 s
   EKF + bias, OCV init : NEES mean 3.48 / per-run median 3.09 (n = 4), runs > 3n: 0,
                          worst run RMS 0.018 %, bias RMS 0.08 mA, median time to
                          stay within 1 % SOC: 0 s
```
