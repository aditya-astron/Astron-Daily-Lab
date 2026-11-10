# 027 - Cabin Humidity Control with a Condensing Heat Exchanger

**Domain:** ECLSS (cabin thermal and humidity control)

## Problem

A crewed module runs a CCAA-style loop: one fan pulls cabin air, a valve
splits it between a condensing heat exchanger (CHX) cooled by a 4.4 degC
water loop and a bypass, and the valve position is set to hold cabin
temperature at 22 degC. Nothing controls humidity directly; it floats to
wherever the coil's condensation balances the crew's vapour output.

For a 0.24 kg/s fan, a CHX with hA = 700 W/K (air side) and 1800 W/K
(coolant side) and 0.075 kg/s coolant:

1. Where do RH and dew point settle for 2 to 8 crew?
2. How does that move with coolant inlet temperature?
3. What is the warmest coolant inlet that still keeps RH <= 60 % with
   6 crew, one exercising continuously (a conservative peak)?

Loads assumed: 1.2 kW equipment + 100 W sensible per crew, 1.9 kg/day
water vapour per crew, +250 W and +0.35 kg/h for exercise.

## Approach

**Psychrometrics.** Hyland-Wexler saturation pressure over water,
`w = 0.621945 pv / (p - pv)`, `h = 1006 t + w (2.501e6 + 1860 t)` J/kg dry air.

**Wet coil (Braun-Klein-Mitchell).** With Lewis number 1, heat and mass
transfer on the air side share one driving potential, `h_air - h_sat(Ts)`.
Linearising `h_sat` with slope `c_s` over the coolant range makes the coil a
dry counterflow exchanger in enthalpy:

    m*   = m_air c_s / (m_cool cp_w)
    NTU* = 1 / (m_air (cp_m / hA_air + c_s / hA_cool))
    Q    = eps_cf(NTU*, m*) m_air (h_in - h_sat(T_cool,in))

`c_s` depends on the coolant outlet, so it is iterated. The outlet air state
follows from an effective saturated surface with `NTU_o = hA_air / (m_air cp_m)`.

**Partly wet coil.** Air meets the warm coolant first, so the coil is dry
over a fraction `phi` of its area and wet over the rest. The boundary is
where the dry-surface temperature reaches the inlet dew point. For each
`phi` the dry (sensible eps-NTU) and wet sections are coupled through the
coolant temperature at the boundary; `phi` is solved from the boundary
condition.

**Cabin balance.** With CHX fraction `f` of the fan flow and a well-mixed
cabin, steady state needs

    f m_fan (h_cabin - h_chx,out) = Q_s + m_v h_g(T)
    f m_fan (w_cabin - w_chx,out) = m_v

Unknowns `f` and `w_cabin`; nested bracketed root solves (Illinois regula
falsi).

**Limits.** Condensate enthalpy is neglected (about 1 % of the latent
duty), the coolant inlet is fixed, the coil has uniform hA, and the
condensate is assumed to be removed without carryover.

## Verification

`test_chx.cpp` (42 checks):

- saturation pressure at 0.01, 20, 25, 40, 60 degC within 0.1 % of IAPWS
  values; 22 degC / 50 % hand calculation of w, h and dew point (11.1 degC)
- dry coil equal to the hand eps-NTU result and to the segmented model
- a 3000-4000 segment counterflow reference coil with no `h_sat` linearisation
  (surface temperature solved locally, coolant outlet found by shooting;
  converged to 1e-5 in Q between 3000 and 12000 segments). Across 36 cases
  (0.03-0.24 kg/s, 3-9 degC coolant, 35-70 % RH) the effectiveness model is
  within 1.6 % on Q and 0.28 K on outlet temperature; at the 4-crew operating
  point within 0.1 % on Q and 0.2 % on condensate
- energy closure: air enthalpy drop = coolant enthalpy rise
- cabin solution: condensate equals the vapour load, CHX duty equals
  `Q_s + m_v h_g`
- 48 h transient of a 100 m^3 cabin with a PI bypass controller from
  25 degC / 70 % RH settles on the algebraic solution (w to 1e-4, valve
  fraction to 1e-3)
- RH rises with crew count; the cabin dew point never falls below the
  coolant inlet

## Run

    make run      # design study
    make test     # verification (about 10 s)

Needs only a C++17 compiler.

## Sample output

```
1) Crew count, coolant inlet 4.4 degC
case          Qs[W]  vap[kg/d]  RH[%]  dew[C]  f_chx   Tout[C]  Tcool,o   Qchx[W]
2 crew          1400     3.80   35.7    6.13   0.332     4.72     9.22     1512
4 crew          1600     7.60   39.0    7.41   0.389     5.14    10.21     1824
6 crew          1800    11.40   42.0    8.50   0.454     5.75    11.20     2135
8 crew          2000    15.20   44.9    9.50   0.531     6.57    12.19     2447
6 + 1 exer      2050    19.80   48.1   10.54   0.564     7.13    12.79     2633

2) Coolant inlet temperature, 4 crew
Tc,in  2.0      1600     7.60   34.2    5.50   0.336     2.45     7.81     1824
Tc,in  8.0      1600     7.60   48.6   10.68   0.520     9.43    13.81     1824
Tc,in 12.0      1600     7.60   64.3   14.96   0.945    15.12    17.81     1824

3) Warmest coolant inlet for RH <= 60 % with 6 crew, 1 exercising
limit           2050    19.80   60.0   13.89   0.900    12.72    16.42     2633
   coolant inlet <= 8.04 degC (dew point 13.89 degC, CHX carries 90 % of fan flow)
```

## Takeaways

- Quadrupling the crew only moves the dew point from 6.1 to 9.5 degC: the
  temperature valve opens to carry the extra sensible load, and that extra
  CHX flow also condenses more water.
- The coolant inlet temperature sets the humidity floor. Each 2 degC
  warmer coolant raises the 4-crew dew point by roughly 1.6-2.2 degC.
- With this hardware the coolant can run up to about 8 degC before the
  exercise peak pushes RH past 60 %, at which point the CHX is already
  taking 90 % of the fan flow, so there is little margin left on either
  temperature or humidity.
