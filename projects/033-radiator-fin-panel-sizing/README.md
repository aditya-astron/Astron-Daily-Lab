# 033 - Minimum-Mass Pumped-Loop Radiator Panel

**Domain:** Thermal analysis

## Problem

A body-mounted radiator has to reject 1.2 kW from a single-phase pumped loop
(50/50 propylene glycol-water, 45 degC in, 15 degC out, 11.1 g/s) to a 220 K
effective sink. The panel is 3/8 in aluminium tubes bonded to aluminium fin
sheet, silvered FEP at end-of-life emissivity 0.80, one radiating face. A
metre of tube with its coolant and bond weighs 0.135 kg. How thick should the
fin sheet be and how far apart should the tubes go to make the panel as light
as possible, and what does a minimum handling gauge or an area budget cost?

## Approach

**Fin.** Straight rectangular fin, root at the tube wall temperature, adiabatic
tip, `n` radiating faces:

    k t T'' = n eps sigma (T^4 - Ts^4)

Multiplying by `T'` and integrating once from the tip (temperature `Tt`) gives
an exact first integral, with `a = n eps sigma / (k t)`:

    T'^2 = 2 a H(T, Tt),   H = (T^5 - Tt^5)/5 - Ts^4 (T - Tt)
    q    = k t sqrt(2 a H(Tb, Tt))
    L    = J / sqrt(2a),   J = int_{Tt}^{Tb} dT / sqrt(H)

`J` becomes a smooth quadrature after `T = Tt + u^2` (which removes the tip
singularity) and `u = sqrt(Tt - Ts) sinh(s)` (which stretches the layer that
forms when a long fin's tip approaches the sink). For a given `L` the tip
temperature comes from Illinois regula falsi on `ln L` against `ln(Tt - Ts)`.
Nothing here is linearised, and the sink is included exactly.

**Minimum-mass fin.** Eliminating `t` and `L` for a required `q` gives the
profile area

    A = t L = q^3 J(Tt) / ((2 n eps sigma)^2 k H(Tb, Tt)^{3/2})

so the best tip temperature minimises `J / H^{3/2}` alone. It does not depend
on `q` or `k`, and the least fin mass scales exactly as `rho / k`.

**Panel.** Each tube carries two fins and its own footprint radiates. With the
wall at the local fluid temperature (film and bond resistances neglected), the
tube length to cool the loop from `Tin` to `Tout` is

    Ltube = int_{Tout}^{Tin} mdot cp / q'(T) dT,   q' = 2 q_fin(T, L) + n eps sigma D (T^4 - Ts^4)

evaluated with 24-point Gauss-Legendre. Panel mass is
`Ltube (m_tube + 2 rho t L)`; Nelder-Mead on `(ln t, ln L)` finds the minimum,
with an optional minimum gauge. For an area budget, the area constraint fixes
`L` at each `t` (area rises monotonically with `L`) and golden section searches `t`.

## Verification (`make test`)

- `J` for a zero sink against the closed form through the incomplete beta
  function, `int_1^th dθ/sqrt(θ^5-1) = [B(3/10,1/2) - B(th^-5; 3/10,1/2)]/5`:
  agreement to 6e-14.
- Tip temperature and root heat flow recovered from the closed form to 1e-11.
- Near-linear fin (`Tb - Ts` = 0.01 K): efficiency matches `tanh(mL)/mL` with
  `m^2 = 4 n eps sigma Ts^3 / (k t)` to 1e-5.
- `k -> infinity` gives efficiency 1.
- An independent finite-difference Newton solver (ghost-node tip) converges
  to the exact heat flow at order 2.000 with the 220 K sink, matches the tip
  temperature to 5e-6 K, and conducts at the root what it radiates.
- Minimum-mass fin, sink 0 K and 220 K: brute force over 301 thicknesses, each
  with the finite-difference length that carries the required `q`, finds the
  same minimum area to 4e-6 and the same thickness to grid resolution. The
  optimum `Tt` is unchanged by `q` and `k`, and area scales as `q^3 / k`.
- Tube length by quadrature in `T` equals an RK4 march of
  `dT/dz = -q'(T)/(mdot cp)` along the tube; marching that length leaves the
  coolant at `Tout` and integrates to exactly 1200 W rejected.
- The panel optimum is stationary (relative gradient 2e-8) and no point on a
  surrounding 21 x 21 grid is lighter.

Hand check of the area bound: at the mean loop temperature 303 K a perfect
fin rejects `0.80 sigma (303^4 - 220^4)` = 276 W/m^2, so 1200 W needs about
4.35 m^2; the integral over the loop gives 4.42 m^2.

## Run

    make run     # design study
    make test    # verification

Needs a C++17 compiler only.

## Sample output

```
Pumped-loop radiator: Q = 1200 W, loop 318.1 -> 288.1 K, mdot = 11.11 g/s
Al 6061 fins, eps = 0.80, one radiating face, sink 220 K
Isothermal (eta = 1) area bound: 4.42 m^2

Designs                        t mm  pitch mm  tube m   area m2  fins kg  total kg
  baseline 1.0 mm, 100 mm      1.00   100.0    45.0    4.50   10.99   17.05
  optimum, no gauge limit      0.11   211.2    34.3    7.24    2.07    6.70
  optimum, t >= 0.5 mm         0.50   212.9    24.4    5.20    6.71   10.01

Mass vs. planform budget (t >= 0.5 mm)
  A_max m2   t mm  pitch mm  tube m  mass kg   eta_fin at inlet/outlet
     4.60    0.50    106.8    43.1    11.47        0.949 / 0.961
     4.80    0.50    149.8    32.0    10.40        0.902 / 0.924
     5.00    0.50    183.4    27.3    10.08        0.860 / 0.890
     5.50    0.50    212.9    24.4    10.01        0.822 / 0.858

Fin material at t >= 0.5 mm (tube line mass unchanged)
  material        k     rho   k/rho    t mm  pitch mm  area m2  mass kg
  Al 1100-O      222    2710  0.0819    0.50    233.4     5.14     9.66
  Al 6061-T6     167    2700  0.0619    0.50    212.9     5.20    10.01
  Mg AZ31B        96    1770  0.0542    0.50    205.5     5.62     8.44
  Cu C10100      391    8940  0.0437    0.50    174.2     4.65    23.25

Minimum-mass single fin (one face, eps 0.80, Tb 300 K, q 50 W/m, Al 6061)
  sink K   Tt/Tb     t mm     L mm    eta
       0  0.7980    0.502    240.7  0.565
     150  0.8182    0.589    253.6  0.572
     220  0.8776    1.123    323.5  0.592
     260  0.9327    3.253    513.2  0.608
```

## Reading the results

- The 1 mm, 100 mm-pitch panel a first sizing might pick is 70 % heavier than
  needed. Tube line mass dominates, so the optimum doubles the pitch to about
  210 mm and accepts fin efficiencies near 0.8.
- Left unconstrained, mass drives the fin to 0.11 mm foil and a 7.2 m^2 panel.
  A 0.5 mm handling gauge costs 3.3 kg but recovers 2 m^2.
- Every area budget is met at the gauge limit by closing the tube pitch, not by
  thickening the sheet. Shrinking from 5.5 to 4.6 m^2 (4 % above the eta = 1
  bound) costs only 1.5 kg.
- `k / rho` ranks fin materials only when thickness is free. With the gauge
  active, density sets fin mass, so magnesium beats both aluminium alloys
  despite its conductivity.
- For a single fin the best tip temperature moves toward the root as the sink
  warms, and the best efficiency stays in a narrow 0.56 to 0.61 band.

## Limits

Steady state, uniform sink, no fin-to-fin or fin-to-tube view factors (flat
panel), no fluid film or bond resistance, and fixed coolant properties. The
first two matter most for a panel that sees a planet or the vehicle body.
