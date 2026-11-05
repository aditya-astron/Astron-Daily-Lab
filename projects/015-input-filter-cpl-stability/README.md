# 015 - Damped LC input filter vs. constant-power-load instability

A 100 W point-of-load converter hangs off a 28 V spacecraft bus behind an LC
filter (10 uH / 22 uF) that is there for conducted-emission limits. A regulated
converter holds its output power constant regardless of input voltage, so its
small-signal input impedance is a negative resistance

    Zin = dV/dI = -V^2 / P          (-7.8 Ohm at 28 V, 100 W)

An undamped LC has a resonant output impedance peak of roughly R0^2/RL
(R0 = sqrt(L/C)), here 22.7 Ohm, far above 7.8 Ohm. Filter plus converter is
then an oscillator and the bus collapses. The goal: size an Rd-Cd damping leg
across the filter capacitor, and verify the fix three independent ways.

## Approach

**Output impedance.** A small complex nodal-analysis solver (stamp R, L, C into
the admittance matrix, inject 1 A at the output, Gaussian elimination) computes
Zout(jw) with the source shorted. It is cross-checked against the hand-derived
expression `1/Zout = 1/(RL + sL) + sC + 1/(Rd + 1/(sCd))`.

**Optimal damping (Erickson, *Fundamentals of Power Electronics*, 10.4.1).**
With n = Cd/C, |Zout| at w = w0 sqrt(2/(2+n)) is the same for Rd = 0 and Rd = inf,
so every Rd curve passes through that point. The optimum Rd puts the peak on it:

    Rd_opt     = R0 sqrt( (2+n)(4+3n) / (2 n^2 (4+n)) )
    |Zout|peak = R0 sqrt( 2(2+n) ) / n

The code also finds Rd by direct numerical minimax (golden section on Rd,
log sweep + golden refinement on frequency) and compares.

**Stability, three ways.**

1. Middlebrook criterion: |Zout| << |Zin| at all frequencies (margin in dB).
2. Routh-Hurwitz on the linearised state equations (iL, vC, vCd) with the
   converter as a conductance -P/V^2. For the bare LC this reduces to
   `RL C > L/Rn`, giving a closed-form power limit used as a test.
3. Nonlinear RK4 simulation of the filter feeding an ideal constant-power load,
   kicked off its operating point by 50-100 mV.

Middlebrook is sufficient but conservative, so the Routh power limit must sit
above the power at which V^2/P equals the peak |Zout|; the test checks that, and
simulates at 0.9x and 1.1x of the Routh limit to confirm the boundary.

## Build and run

    make
    ./test_filter     # verification suite
    ./filter_demo     # design report

Needs a C++17 compiler; no dependencies.

## Sample output

```
Bus 28 V, converter 100 W, L = 10 uH, C = 22 uF, RL = 20 mOhm
R0 = sqrt(L/C) = 0.674 Ohm, f0 = 10.73 kHz
Converter input: v = 27.928 V, Zin = -V^2/P = -7.80 Ohm

Damping leg (Erickson optimum, n = Cd/C = 4): Cd = 88 uF, Rd = 0.4129 Ohm
Numerical minimax Rd (with RL included):         Rd = 0.3934 Ohm

filter     peak|Zo| Ohm   f_peak kHz   Middlebrook dB    Routh        P_max W     sim growth x
undamped         22.737        10.73             -9.3 UNSTABLE           34.4     bus collapse
damped            0.546         6.37             23.1   stable         1384.8         7.03e-12
```

The bare filter can only feed 34 W before it goes unstable. An 88 uF / 0.41 Ohm
leg (an electrolytic with that ESR often does the job by itself) drops the peak
to 0.55 Ohm, gives 23 dB of Middlebrook margin at 100 W and pushes the stability
limit past 1.3 kW.

## Verification

- Nodal solver vs closed-form Zout over 100 Hz-1 MHz: max relative error ~1e-15.
- Bare LC resonance at f0 = 1/(2 pi sqrt(LC)) and peak R0^2/RL to 0.1%.
- Numerical minimax Rd, peak |Zout| and peak frequency match Erickson's
  closed forms for n = 1, 4, 10.
- Undamped Routh limit matches the closed form P = v^2 RL C / L,
  v = Vin / (1 + RL^2 C / L); simulation decays at 0.8x and grows at 1.2x.
- Damped design: Routh and simulation agree at 100 W and on both sides of
  the power limit.

## Limits

The converter is an ideal CPL at all frequencies. A real one only looks like
-V^2/P below its control-loop crossover, so this is the worst case; above
crossover its input impedance turns positive and margin improves. Capacitor
ESR and ESL of C are neglected (ESR in C helps damping).
