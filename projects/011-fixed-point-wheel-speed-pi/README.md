# 011 - Fixed-point PI speed loop for a reaction wheel with anti-windup

## Problem

A small-satellite reaction wheel runs a speed loop on a low-cost MCU with no
FPU. The goal is an integer-only PI controller (Q16.16) that

- gives a designed closed-loop response (wn = 8 Hz, zeta = 0.8),
- survives large speed commands where the 1.5 A motor driver saturates,
- holds zero steady-state error under a load torque,

and to show that the integer implementation actually delivers the response
the continuous-time design promises.

## Approach

**Plant.** With a fast inner current loop, `J dw/dt = Kt i - b w - tau`.
Current and load are held for each 1 ms sample, so the plant is advanced with
its exact solution `w = w_inf + (w0 - w_inf) exp(-b dt / J)`, with no
integrator error in the reference.

**Gain design.** With `i = Kp e + Ki int(e)` the closed loop is
`J s^2 + (b + Kt Kp) s + Kt Ki`. Matching `J (s^2 + 2 zeta wn s + wn^2)`:

    Kp = (2 zeta wn J - b) / Kt      Ki = wn^2 J / Kt

**Analytic references.**

- Speed step: `W/R = (a s + wn^2) / (s^2 + 2 zeta wn s + wn^2)`, `a = Kt Kp / J`,
  i.e. the standard second-order step plus `a` times its impulse response.
- Current-limited slew: `t = (J/b) ln(w_inf / (w_inf - w1))`, `w_inf = Kt Imax / b`.
- Load step: `W/Td = -s / (J (s^2 + 2 zeta wn s + wn^2))`, so the dip is
  `-(tau / (J wd)) e^{-sigma t} sin(wd t)`, peaking at `wd t = atan(wd / sigma)`.

**Fixed point.** Q16.16 in int32, products in int64 with round-to-nearest
(truncation biases the integrator by up to 1 LSB per sample). Anti-windup is
conditional integration: the integrator freezes while the output is at the
limit and the error pushes further into it, and it is bounded to +/-Imax.

One consequence of fixed point worth knowing: the per-sample integrator
increment `Ki Ts e` rounds to zero for `|e| < 0.5 / (Ki Ts)` = 19.8 LSB, so
the loop has a dead band of about 3e-4 rad/s around the setpoint. That is the
residual error seen under load; widening the integrator accumulator (e.g. a
Q8.24 state) would shrink it if it mattered.

## Files

| File | Purpose |
|------|---------|
| `pi_q16.c/.h` | Integer PI controller, the part that goes on the MCU |
| `wheel.c/.h` | Exact-discretized wheel dynamics |
| `loop.c/.h` | Gain design, analytic responses, closed-loop harness, double-precision twin |
| `wheel_demo.c` | Small step, saturated step with/without anti-windup, load step |
| `test_loop.c` | 15 checks against hand calculations and analytic responses |

## How to run

    make            # builds wheel_demo and test_loop (gcc, libm)
    make test
    ./wheel_demo

## Sample output

```
Design: wn=50.27 rad/s zeta=0.80 -> Kp=0.8042 A/(rad/s)  Ki=25.266 A/rad
Q16.16: Kp=52704  Ki*Ts=1656 (LSB)

1) 1 rad/s step at 100 rad/s, sampled loop vs analytic continuous
   max |sim - analytic| = 0.0188 of step

2) 0 -> 300 rad/s (current-limited, analytic t90 = 1.808 s)
   no AW    overshoot  97.87 %   t90 1.809 s   not settled in 10 s (saturated limit cycle)
   clamp AW overshoot   0.11 %   t90 1.809 s   1% settle 1.990 s

3) 5.0 mN m load step at 300 rad/s
   peak dip -0.2168 rad/s (analytic -0.2109), final error -6.75e-05 rad/s
```

The small-signal mismatch (1.9 % of the step) and the 3 % deeper load dip are
the expected cost of sampling at 1 kHz against an 8 Hz loop, not quantization:
the integer controller tracks its double-precision twin to within 8 LSB.

Without anti-windup the wheel does not merely overshoot. Because the plant is
almost a pure integrator (`J/b` = 200 s), the wound-up integrator drives a
sustained relay-like oscillation between roughly 20 and 600 rad/s that never
dies out. Clamping the integrator removes it entirely.
