# 044 - Attitude control of a spacecraft with flexible solar arrays

## Problem

A 1740 kg m^2 spacecraft (600 kg m^2 hub, two 6 m, 30 kg panels) runs a PD
attitude loop at 2 Hz with a 0.2 rad/s rigid-body bandwidth. The panels'
first antisymmetric bending mode sits at 0.40 Hz, only 12x above the
bandwidth and 1/5 of the sample rate. Does the loop survive it, and if not,
which fix is robust to the 10-30 % error typical of predicted panel
frequencies? Candidates: a notch on the mode, a notch below the mode, a
low-pass roll-off (gain stabilisation), and simply sampling at 10 Hz.

## Approach

**Structure.** Each panel is an Euler-Bernoulli cantilever on a hub of radius
r, discretised with 40 Hermite cubic elements. Panel points move as
y = (r + x) theta + w(x), with w clamped at the root, which gives the hybrid
(constrained-mode) equations

    J_t theta'' + sum L_k eta_k''                     = tau
    L_k theta'' + eta_k'' + 2 z w_k eta_k' + w_k^2 eta_k = 0,   L_k = 2 int m (r + x) phi_k dx

For one mode the hub transfer theta/tau has zeros at the cantilever frequency
w_c and poles at the free-free frequency w_p = w_c / sqrt(1 - L^2/J_t);
sensor and actuator are collocated on the hub, so poles and zeros alternate.
Three modes carry 98.8 % of the panel inertia (sum L_k^2 -> panel inertia).
Structural damping 0.5 %.

**Loop.** Hub angle and gyro rate sampled at t_k, torque held (ZOH) from
t_{k+1}: one sample of transport delay. Kp = J_t wn^2, Kd = 2 zeta wn J_t.
Filters are biquads by prewarped bilinear transform. The plant is
discretised exactly (matrix exponential) and closed with the filter and delay
states into one augmented matrix; stability is its spectral radius. Margins
come from the discrete Nyquist plot, frequency tolerance from bisecting the
unit-circle crossing while scaling every modal frequency.

**Why it goes unstable.** Rate feedback on a collocated mode adds damping only
while the PD lead at w_p, atan(2 zeta w_p / wn), exceeds the sampling lag
1.5 w_p T (half sample ZOH + one sample delay). At w_p = 2.51 rad/s that
allows T < 0.40 s; 2 Hz is past it. A notch lags below its centre and leads
above it, so it only phase-stabilises a mode sitting at or above the notch.

## Verification

`python3 -m pytest -q test_flex_attitude.py` (16 tests, ~17 s):

- FE cantilever frequencies vs beta^2 sqrt(EI/(m l^4)), beta = 1.8751, 4.6941 (1e-7, 1e-6).
- L_1 vs quadrature of the exact clamped-free mode shape (1e-6); sum L_k^2
  equals the panel inertia (Parseval) to 1e-3.
- Hybrid model with all modes reproduces the full hub+panel FE free-free
  frequencies to 1e-9; one-mode model equals the closed-form w_p.
- ZOH discretisation and the full closed loop (filter, delay, ordering) vs a
  sample-by-sample continuous solve_ivp run to 1e-10.
- Rigid loop at T = 1 ms recovers the continuous wn and zeta to 0.2 %.
- Prewarped notch is exactly zero at w0; low-pass is exactly 1/(2 zeta), -90 deg at wf.
- Gain margins from Nyquist phase crossings equal the eigenvalue-sweep gain
  limits to 0.2 % for three designs.
- Zero steady-state error and zero holding torque from the closed-loop fixed point.
- Frequency-tolerance edge unchanged going from 3 to 5 modes.
- 1.5 w T phase rule within 5 % of the exact critical period at nominal.

## Run

    pip install numpy scipy pytest
    python3 flex_attitude.py
    python3 -m pytest -q test_flex_attitude.py

## Output

```
Structure
  J_t = 1740 kg m^2 (hub 600, panels 1140)
  mode 1: cantilever   1.575 rad/s, L =  32.50 kg^0.5 m, free-free   2.508 rad/s (0.399 Hz)
  mode 2: cantilever   9.869 rad/s, L =   7.58 kg^0.5 m, free-free  10.329 rad/s (1.644 Hz)
  mode 3: cantilever  27.634 rad/s, L =   3.48 kg^0.5 m, free-free  27.909 rad/s (4.442 Hz)
  modal effective inertia (panels): 1 mode 0.9266, 2 0.9770, 3 0.9876
  one-mode estimate w_c / sqrt(1 - L1^2/J_t) = 2.5124 rad/s vs FE 2.5079
  at 2 Hz, modes 2 and 3 fold to 2.24 and 2.78 rad/s

Rigid-body PD: wn = 0.2 rad/s, zeta = 0.7, Kp = 69.6 N m/rad, Kd = 487.2 N m s/rad

design                           stable   worst flex pole rigid PM      GM  w_p range (x nominal)
-------------------------------------------------------------------------------------------------
PD, 2 Hz                             NO -0.0348 @ 2.69 r/s    51.9d       -                      -
PD + notch at w_p, 2 Hz             yes +0.0082 @ 2.50 r/s    43.4d  13.9dB            0.95 - 1.36
PD + notch at 0.8 w_p, 2 Hz         yes +0.0240 @ 2.24 r/s    40.7d  10.3dB            0.76 - 1.24
PD + low-pass 0.8 rad/s, 2 Hz       yes +0.0085 @ 2.49 r/s    27.1d   4.9dB            0.77 - 2.50
PD, 10 Hz                           yes +0.0050 @27.91 r/s    62.1d  16.9dB            0.30 - 2.50
  worst flex pole: damping ratio and frequency of the least-damped pole above 1 rad/s (open loop 0.005)
  GM: loop-gain increase to the first unit-circle crossing (gain reductions stay stable down to 0.001x)

Slew response, 0.5 deg step
  PD, 2 Hz                         diverges: tip 260.9 mm after 60 s
  PD + notch at w_p, 2 Hz          settle(2%)   25.0 s, peak torque 0.589 N m, peak tip 1.20 mm, last 20 s 0.001 mm
  PD + notch at 0.8 w_p, 2 Hz      settle(2%)   24.5 s, peak torque 0.553 N m, peak tip 1.26 mm, last 20 s 0.000 mm
  PD + low-pass 0.8 rad/s, 2 Hz    settle(2%)   18.5 s, peak torque 0.609 N m, peak tip 1.17 mm, last 20 s 0.006 mm
  PD, 10 Hz                        settle(2%)   29.4 s, peak torque 0.607 N m, peak tip 1.82 mm, last 20 s 0.000 mm

Longest stable sample period for plain PD vs the 1.5 w T phase rule
  w_p x 0.6: exact 0.5727 s, phase rule 0.6540 s (+14.2 %)
  w_p x 0.8: exact 0.4606 s, phase rule 0.4983 s (+8.2 %)
  w_p x 1.0: exact 0.3857 s, phase rule 0.4024 s (+4.3 %)
  w_p x 1.2: exact 0.3325 s, phase rule 0.3375 s (+1.5 %)
  w_p x 1.5: exact 0.2766 s, phase rule 0.2716 s (-1.8 %)
```

## Takeaways

- Plain PD at 2 Hz is unstable: the first panel mode goes to -3.5 % damping.
  Its critical sample period is 0.386 s; the phase rule gives 0.402 s.
- A notch centred on the predicted mode works nominally but fails if the real
  mode is more than 5 % low, where the notch's phase lag adds to the delay.
- Moving the same notch to 0.8 w_p keeps the mode in the notch's lead region:
  stable from 0.76 to 1.24 x nominal with 41 deg rigid phase margin. Above
  1.24 x the mode leaves the notch and the delay wins again.
- The 0.8 rad/s low-pass is stable for any stiffer-than-0.77 x mode but costs
  rigid phase margin (27 deg, 4.9 dB GM): an underdamped slew.
- At 10 Hz plain PD actively damps the mode (9.5 % damping on mode 1) and
  tolerates any mode frequency from 0.3 to 2.5 x: when the computer allows
  it, sampling faster beats filtering.
- At 2 Hz modes 2 and 3 fold to 2.24 and 2.78 rad/s, right beside mode 1.
  They are kept in the model; any filter designed on the continuous picture
  alone would miss them.
