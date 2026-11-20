# 050 - Regenerative cooling of a LOX/methane thrust chamber

A 7 MPa LOX/methane chamber with a 60 mm throat (about 10.7 kg/s, 2.4 kg/s
of fuel) is cooled by running all of the methane through milled channels in
a CuCrZr liner, counterflow from the nozzle exit (area ratio 8) to the
injector. The questions: where the hot wall peaks, whether a 900 K hot-wall
limit and a 25 % of Pc pressure-drop budget can both be met, and which
channel count and height taper does it with the least pressure drop.

## Approach

**Contour.** Cylinder, fillet, 30 deg cone, throat arcs of 1.5 Rt upstream
and 0.382 Rt downstream, 15 deg cone. Contraction ratio 4; the cylinder
length comes from L* = 1.0 m.

**Gas side.** Frozen products, T0 = 3550 K, gamma = 1.15, M = 21 kg/kmol
(inputs of the order an equilibrium code gives at O/F 3.4). Local Mach from
the isentropic area ratio, turbulent recovery temperature with r = Pr^(1/3),
and Bartz:

```
hg = 0.026/Dt^0.2 (mu^0.2 cp / Pr^0.6) (pc/c*)^0.8 (Dt/rc)^0.1 (At/A)^0.9 sigma
sigma = 1 / [ (Twg/(2 T0) (1 + (g-1)/2 M^2) + 1/2)^0.68 (1 + (g-1)/2 M^2)^0.12 ]
```

with Bartz's own mu and Pr estimates. sigma depends on Twg, so each station
is an implicit solve.

**Wall and channel.** Per channel, gas heat over the hot arc 2 pi r / N goes
through the ligament tw and leaves through the channel floor and two land
faces, treated as adiabatic-tip fins of thickness b:

```
q (2 pi r / N) = hc (Twc - Tb)(w + 2 eta h),  eta = tanh(mh)/(mh),  m = sqrt(2 hc/(k b))
Twc = Twg - q tw / k
```

Coolant hc is Gnielinski with a Sieder-Tate (mu_b/mu_w)^0.14 correction,
friction is Haaland with 1.6 um roughness, and methane properties come from
CoolProp (HEOS) on (h, p). Land width is fixed at 0.8 mm; channel height
tapers linearly in radius from h_t at the throat to h_c at the chamber
radius. Brent finds Twg at each station; a Heun march carries (h, p) along
the flow with ds/dx = sqrt(1 + r'^2) on heat input and friction length.

**Verification** (`test_regen.py`, 19 tests)

- Area-Mach inversion against NACA 1135 (A/A* = 2: M = 0.3059 and 2.1972),
  round trips at gamma = 1.15, c* against the hand formula.
- Bartz in SI equals the original Huzel & Huang inch-pound-Btu form to 2e-4.
- Contour value and slope continuous at every joint; integrated volume
  reproduces L* to 1e-5.
- Constant-property cylinder: coolant outlet temperature matches the
  exponential `Taw - Tb(L) = (Taw - Tb0) exp(-U'L / mdot cp)` to 1 mK, the
  pressure drop matches `f L/Dh G^2/(2 rho)`, and the march error falls 4x
  per halving of step (second order).
- Real chamber: total heat equals an independent trapezoid of the recorded
  flux to 0.2 %, outlet temperature converges under refinement, station
  balance closes to 1e-8.
- 2D finite-volume conduction of one channel/land cell: exact in the slab
  limit, mesh-convergent, and within 2 % of the fin model's heat.

## Run

```
pip install numpy scipy CoolProp pytest
python3 run_regen.py
python3 -m pytest -q test_regen.py
```

## Results

```
Baseline jacket: N = 120, h = 3.0 mm at throat -> 1.5 mm at chamber, land 0.8 mm, ligament 0.7 mm
     x[mm]  r[mm]     M   hg[kW/m2K]  q[MW/m2]  Twg[K]  Twc[K]  Tb[K]  p[MPa]
     206.2   84.9  3.00       3.65     11.23   399.3   374.7  120.0  12.000
       0.5   30.0  1.03      23.95     64.05   863.0   722.9  292.6  11.791
     -20.2   34.8  0.51      18.60     49.64   877.9   769.3  325.8  11.673
     -60.3   57.6  0.16       7.65     21.27   771.5   725.0  385.2  11.550
    -280.7   60.0  0.15       6.99     18.63   885.2   844.5  600.0  10.676
  peak q 64.10 MW/m2 at x = -0.7 mm
  heat load 4.034 MW, coolant 120 -> 600.0 K, dp 1.324 MPa

2D unit cell at the throat vs. the 1D fin model on the same flat cell
  1D: Twg 871.2 K
    n   Q2D/Q1D   Twg max 2D   Twg min 2D
   80    0.9954        887.1        879.9

  lowest-dp feasible: N = 132, h = 3.0/1.8 mm, dp 1.161 MPa, Twg max 896.8 K
```

What the numbers say:

- Heat flux peaks at 64 MW/m2 just upstream of the throat, but the hottest
  wall is at the injector end. By then the methane has absorbed 4 MW and is
  at 600 K, so the film drop has little room. Most of the trade's feasible
  designs are limited there, not at the throat.
- Deeper throat channels lower the pressure drop but heat the throat; shallow
  chamber channels raise the coolant velocity where the methane is hot and
  gas-like. 96 channels fail at the throat whatever the taper; 132 channels
  at 3.0/1.8 mm meet both limits with the least pressure drop (1.16 MPa).
- The 1D fin model gets the heat right to 0.5 %, but the 2D hot face runs
  9 to 16 K above it because the ligament over the channel centre is fed
  only sideways. That is the margin to hold against the 900 K limit.

Limits: Bartz with a fully developed boundary layer overpredicts heat flux
in the barrel, which makes the injector-end result conservative. Heat
transfer deterioration near the pseudo-critical point, channel curvature and
thermal stratification are not modelled. No film cooling.
