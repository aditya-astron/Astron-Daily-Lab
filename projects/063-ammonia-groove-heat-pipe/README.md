# 063 - Axially Grooved Aluminium/Ammonia Heat Pipe: Transport Limits and Groove Sizing

**Domain:** Thermal analysis

## Problem

An 80 W avionics box is tied to a radiator by a 12.7 mm OD aluminium/ammonia
heat pipe with axial rectangular grooves: evaporator 0.15 m, adiabatic 0.60 m,
condenser 0.30 m, vapour core (fin-tip circle) 8.5 mm, fins at least 0.35 mm
wide, wall under the grooves at least 0.5 mm (so d <= 1.6 mm). It must carry
80 W with 1.5x margin in orbit from -20 to +60 C, and it must be acceptance
testable on the ground at the usual 2.54 mm adverse tilt with a dryout power
of at least 40 W. What groove width and depth, and how far can a ground test
be trusted for the 0-g number?

Python (numpy, scipy, CoolProp for saturated ammonia).

## Approach

**Capillary limit.** Flat meniscus at the condenser end cap, tightest one
(radius w/2, cylindrical) at the evaporator end cap:

```
2 sigma / w >= Q Leff (R_l + R_v) + rho_l g (h_tilt + Dv)       (Dv term only in 1 g)
R_l = mu_l / (rho_l hfg N C(w, d)),     R_v = (fRe/16) 128 mu_v / (pi rho_v hfg Dv^4)
Leff = La + (Le + Lc)/2   (uniform heating)
```

**Groove flow.** An open groove with a shear-free surface is half of a
symmetric w x 2d duct, so `C = Q mu / G` comes from the exact rectangle series
(White, eq. 3-48). It is checked against a Richardson-extrapolated
finite-difference Poisson solve on the groove itself (Neumann at the free
surface), the Shah & London table and fit, and the deep-slot limit
`w^3 d / 12 (1 - 0.630 w / 2d)`.

**Vapour.** Hagen-Poiseuille below Re 2300, Blasius above. Because the drop
jumps at transition, the limit is found as the largest Q with
`Q (R_l + R_v(Q)) Leff <= drive` rather than by naive iteration (which can
cycle at the boundary).

**Other limits** (Faghri ch. 4): Busse sonic and viscous, Weber entrainment
with wick hydraulic radius w/2. Boiling is not modelled; axial grooves with a
thin evaporator wall run far below it at these fluxes.

**Effective length.** An ODE march of the axial heat flow and meniscus
pressure for any heat distribution, so `Leff = (1/Q) int Qx dx` is checked for
uniform heating and for a box mounted on either end of the saddle.

**Ground test.** In 1 g the top grooves need an extra `rho_l g Dv`, which is
the same as Dv of extra tilt. Dryout power is linear in tilt while laminar, so
a fit over a few small tilts evaluated at tilt = -Dv is the 0-g claim. The
noise gain of that extrapolation is `sqrt(x0' (X'X)^-1 x0)`.

## Run

```
pip install numpy scipy CoolProp pytest
python study.py
pytest -q
```

## Results (sample_output.txt)

- Series and FD groove conductance agree to 1e-7; exact fRe reproduces the
  Shah & London table to 1e-3 (their polynomial fit is off by up to 0.06 %).
- Hand calculation with NIST table properties: 185.1 W for the 0.5 x 1.2 mm
  baseline at 20 C; the model with CoolProp properties gives 185.7 W.
- Baseline: 123 W at 60 C, so it misses 1.5x. Ammonia's transport factor falls
  from 1.31e11 W/m^2 at -20 C to 7.6e10 at 60 C, so the hot end sizes the
  pipe. Sonic and viscous limits are 2-6 orders away; entrainment is 3-10x.
- The vapour core runs at Re ~2400-2600 at the 0-g limit from -40 to 20 C:
  right on transition. It is only a few % of the drop but decides whether
  -20 C or +60 C is the worst case for the wider grooves.
- Width trade at d = 1.6 mm: 0-g capacity rises with w (flow ~ N w^2 d for
  deep grooves) up to ~1.2 mm, but the 1-g static height `2 sigma/(rho g w) - Dv`
  goes to zero at w ~ 0.85 mm. Wide grooves are untestable on the ground.
- Selected: w = 0.55 mm, d = 1.6 mm, N = 29. Worst-case 0-g 188 W (2.35x, at
  60 C, capillary-limited), ground test 46 W at 2.54 mm.
- Box placement alone moves Leff from 0.765 to 0.885 m (+-7 % capacity).
- Ground-test extrapolation from 0-3.8 mm tilts over-predicts the 0-g limit
  by 3.7 %: tests run at Re ~1300, the 0-g limit at Re ~3600 (turbulent core).
  The lever arm to -8.5 mm gives a noise gain of 3.7, so +-3 W per dryout
  point becomes +-11 W on the flight claim.

```
5. Design selection
-------------------
  11 of the 26 x 13 grid meet 0-g >= 120 W (-20..60 C) and 1-g >= 40 W at 2.54 mm
  selected w = 0.55 mm, d = 1.6 mm, N = 29: worst 0-g 188 W (2.35x), ground 46 W
  worst case at 60 C, binding limit: capillary_0g

6. Ground test: Q_max vs adverse tilt at 20 C, extrapolated to 0 g
------------------------------------------------------------------
  fit -21.89 W/mm, zero-flow tilt 4.64 mm; at tilt = -Dv: 287.7 W (model 0-g at 20 C: 277.3 W)
  the +3.7 % bias is the vapour core: Re 1324 in the tests, 3612 (turbulent) at the 0-g limit
  noise gain at -Dv: 3.70 (vs 0.50 inside the data); +/-3 W per dryout point -> +/-11 W on the 0-g claim
```

## Limits of the model

Flat-meniscus liquid area (the receding meniscus near the evaporator trims
some flow area), no vapour-liquid shear on the groove surface, no
liquid-puddle flow in 1 g (makes the ground numbers conservative), and
isothermal properties along the pipe.
