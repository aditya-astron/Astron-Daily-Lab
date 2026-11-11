# 030 - Buck Converter: Exact Periodic Steady State, DCM Boundary and Peak-Current-Mode Stability

**Domain:** Electronics and circuit analysis

## Problem

A point-of-load buck converter takes the 28 V spacecraft bus down to 5 V at
5 A (500 kHz, 10 uH, 47 uF). Before committing to a layout, the averaged
design equations need to be checked against the real switched circuit:
output voltage with conduction losses, inductor and output ripple, efficiency,
and where the converter drops into discontinuous conduction at light load.
A second, 5 V to 3.3 V rail (D = 0.66) runs peak current mode, which is
unstable above 50 % duty unless a compensating ramp is added; find the ramp
that fixes it and show what happens without it.

## Approach

**Exact switched model.** With `x = [iL, vC]` the circuit is linear inside
each switch state, `dx/dt = A x + B u`:

    L diL/dt = u Vin - (ron + rL) iL - vo
    C dvC/dt = (R iL - vC) / (R + rc),   vo = R (vC + rc iL) / (R + rc)

Each interval has the exact solution `x(t+tau) = Phi x + Gamma u`, both taken
from one augmented matrix exponential `expm([[A, B u],[0, 0]] tau)`. The
one-period map is `x -> Phi_off Phi_on x + Gamma`, and its fixed point
`(I - Phi) x* = Gamma` is the periodic steady state, solved in one linear
solve. Period averages use Gauss-Legendre quadrature on each interval.

**Averaged references (synchronous CCM).**

    Vo  = D Vin R / (R + ron + rL)
    dIL = (Vin - Vo - (ron + rL) Io) D T / L
    dVo = dIL / (8 fs C)  (capacitor)   and   dIL rc  (ESR)
    eta = R / (R + ron + rL)

**DCM (non-synchronous, ideal diode).** The diode interval ends when `iL`
reaches zero; that instant is found by root-finding on the exact solution,
then the idle interval holds `iL = 0`. The PSS is the root of the nonlinear
period map. Reference: `M = 2 / (1 + sqrt(1 + 4K/D^2))`, `K = 2L/(R T)`,
with DCM for `K < 1 - D`.

**Peak current mode.** With a stiff output the valley current obeys a 1-D
map. The switch turns off when `iL = Ic - ma t`; with up-slope
`m1 = (Vin - Vo)/L` and down-slope `m2 = Vo/L`, a perturbation scales by

    alpha = -(m2 - ma) / (m1 + ma)

every cycle. Stability needs `ma > (m2 - m1)/2`; `ma = m2` is deadbeat.
Ridley's sampled-data model puts a double pole at fs/2 with
`Q = 1 / (pi (mc D' - 0.5))`, `mc = 1 + ma/m1`, which diverges at the same
ramp. Without a ramp the map is piecewise linear (slope `alpha` when the
comparator trips, slope 1 when the duty clamps), so its Lyapunov exponent is
the unclamped fraction times `ln|alpha|`.

**Limits.** Ideal switching (no dead time, gate charge or diode drop),
linear L and C, stiff output in the current-mode study (no outer voltage loop).

## How to run

    pip install numpy scipy pytest
    python3 buck.py
    python3 -m pytest -q test_buck.py

## Verification

`test_buck.py` checks:

- matrix-exponential PSS equals 600 periods of brute-force DOP853 integration
  from rest (rtol 1e-6)
- Vo = 5.000 V at the trimmed duty, inductor ripple within 0.1 % of the
  formula, output ripple between the ESR and cap terms and their sum
- with rc = 0 the output ripple matches `dIL/(8 fs C)` within 0.1 %
- energy balance: `Pin - Pout` equals the summed `I^2 R` losses to 1e-9,
  and inductor RMS matches `Io^2 + dIL^2/12`
- DCM conversion ratio against the closed form for three loads (2e-4), and
  the mode flips across `R_crit = 2 L fs / (1 - D)` (+/-3 %)
- current-mode perturbation ratio equals `alpha` for four ramps; no ramp
  grows by `|alpha|` per cycle then stays bounded with positive Lyapunov
  exponent; `ma = m2/2` converges; `|alpha| = 1` and Ridley `Q -> inf` give the
  same ramp; deadbeat gives `Q = 2/pi`

Hand check: `dIL = 22.825 V x 0.18482 x 2 us / 10 uH = 0.8437 A`,
`dVo = 0.8437 / (8 x 500 kHz x 47 uF) = 4.488 mV`.

## Sample output

```
1) Synchronous buck 28 V -> 5 V, 5 A, 500 kHz, exact PSS
   D = 0.18482
   Vo avg     PSS 5.00000 V   averaged 5.00000 V
   dIL p-p    PSS 0.84361 A   analytic 0.84371 A
   dVo p-p    PSS 9.003 mV  cap 4.488 mV + ESR 8.437 mV
   efficiency PSS 96.6084 %   averaged 96.6184 %
   energy balance  Pin - Pout = 877.659 mW,  sum I^2 R = 877.659 mW
   rc = 0:  dVo p-p PSS 4.4886 mV   dI/(8 fs C) 4.4878 mV

2) Non-synchronous lossless buck, light load, D = 0.25
   K_crit = 1 - D -> R_crit = 13.333 ohm
   R =    5.0  mode CCM  M sim 0.25000  analytic 0.25000
   R =   20.0  mode DCM  M sim 0.29655  analytic 0.29654
   R =   50.0  mode DCM  M sim 0.42422  analytic 0.42419
   R =  200.0  mode DCM  M sim 0.65589  analytic 0.65587

3) Peak current mode, 5 V -> 3.3 V POL (D = 0.66), 1 MHz
   no slope comp              alpha -1.9412 (map -1.9412)  Ridley Q -1.99
   ma = (m2-m1)/2 boundary    alpha -1.0000 (map -1.0000)  Ridley Q +inf
   ma = m2/2                  alpha -0.4925 (map -0.4925)  Ridley Q +1.87
   ma = m2 deadbeat           alpha -0.0000 (map -0.0000)  Ridley Q +0.637
   no slope comp, cycles 195-200 valley currents: 1.663 2.436 2.594 2.287 2.884 1.726 A
   Lyapunov exponent +0.445 per cycle (> 0: bounded, aperiodic, chaotic)
   with ma = m2/2: -0.708 per cycle
```
