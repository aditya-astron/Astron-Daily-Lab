# 066 - Hot-Swap Inrush Control and MOSFET Junction Temperature

**Domain:** Electronics and circuit analysis / Spacecraft power systems

## Problem

A 120 W payload converter with a 470 uF input bank is switched onto a 28 V
spacecraft bus through a high-side N-MOSFET whose output dV/dt is set by a
gate capacitor (10 uA charge pump, gate clamp 12 V above the bus). The FET
sits on a 60 C cold plate and its junction is derated to 110 C. Which gate
capacitor, and should the converter be allowed to start at its 18 V UVLO
during the ramp or be held off until power-good?

Python (numpy, scipy).

## Approach

**Circuit.** Gate and output nodal equations with a level-1 square-law FET
(Vth 3 V, k 20 A/V^2, Cgs 3 nF):

```
Ig           = Cg dVg + Cgs (dVg - dVo)
Id(Vgs, Vds) = C_L dVo + Vo/R_L + [P/Vo] - Cgs (dVg - dVo)
```

integrated with Radau in three phases (converter off, converter on, gate
clamped) split by events so the load step is resolved exactly.

**Ramp closed forms.** With the output slewing at S = Ig/Cg from 0 to Vb:

```
E_C = C_L Vb^2 / 2                          exact for any profile and any S
E_R = Vb^3 / (6 R_L S)
E_P = (P/S) [Vb ln(Vb/V_en) - (Vb - V_en)]  converter on from V_en during the ramp
```

**Turn-on stall.** When the converter starts, the FET current jumps by P/V_en
and Vgs has to rise by `sqrt(2 I+/k) - sqrt(2 I-/k)` while the gate keeps
slewing at S. The output stalls for about `dVgs/S` with `(Vb - V_en)` across
the FET, which adds `(Vb - V_en)(P/V_en) dVgs / S`.

**Junction temperature.** Five-term Foster network (Zth 4.0 K/W junction to
plate). For power linear between samples each Foster state has an exact
one-step update; it is written with `expm1` so near-zero steps (a power step)
do not cancel.

**Gate timing.** Power-good is the gate reaching its clamp:
`t = Vb Cg/Ig + Vgs_clamp (Cg + Cgs)/Ig`, since Cgs loads the gate whenever
the output is still.

## Run

```
pip install numpy scipy pytest
python study.py
pytest -q
```

## Results (sample_output.txt)

- Capacitive-only charge: simulation matches `C Vb^2/2 = 0.18424 J` to 2e-6
  at 100 V/s through 20 kV/s. The resistive and power-good-gated cases match
  the ramp formulas within 0.4 %.
- Starting the converter at UVLO costs 8-15 % more than the ideal ramp
  formula. The Vgs step estimate explains 90-96 % of that excess at slow
  slews (100-300 V/s). At fast slews the C_L dip and recharge after turn-on dominate. The
  dip stays under 0.8 V, inside a normal UVLO hysteresis, so the latched
  converter model holds.
- Foster stepping reproduces `P Zth(t)` to 1e-9 and agrees with a Radau solve
  of the network on the real ramp waveform to 1e-7 K.
- UVLO start has an interior optimum: slow ramps keep the FET in its linear
  region under converter current (E_P ~ 1/S), and fast ramps push the
  capacitive power `C S Vb` into a short pulse. For 120 W the best slew is
  ~2.1 V/ms with Tj 83.6 C. The 110 C limit is reached near 330 W, and only
  at ~19 V/ms with a 26 A inrush, so in practice the inrush limit binds well
  before that.
- With power-good gating, the ramp heating does not depend on converter power
  (1 kW gives the same 66.3 C at 1 V/ms as 120 W does). Tj falls
  monotonically with slower slew, so the slowest ramp that fits the start-up
  budget wins.
- Pick: **Cg = 15 nF**, power-good gated. With Ig +-20 % and Cg +-10 %, power-good
  comes at 87 ms at the slow corner (22 nF misses the 100 ms budget). At the
  fast corner (0.89 V/ms), Tj is 65.8 C, the inrush peak is 4.32 A (the load
  itself) and the FET absorbs 0.19 J. The same part at the same corner
  with UVLO start reaches 86.0 C, 7.33 A and 0.57 J.

## Limits

Square-law FET with constant Cgs (no Miller nonlinearity in Cgd, no
temperature dependence of Vth or k), ideal bus source, converter modelled as
an ideal constant-power sink from turn-on. Foster parameters are
representative of a D2PAK on a board-to-plate path, not a specific part.
