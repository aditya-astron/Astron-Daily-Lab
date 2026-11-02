# 007 - Transient Orbital Thermal Model of a Nadir-Pointing 3U CubeSat

**Domain:** Thermal analysis

## Problem

A 3U CubeSat flies nadir-pointing in a 500 km circular orbit with 4 W of
average electronics dissipation. Body-mounted cells cover the +/-Y and zenith
faces; the rest is anodised aluminium. How hot and cold do the panels and the
internal stack get over an orbit, and how does that change with the solar
beta angle, from a 38% eclipse at beta = 0 to full sun above about 68 deg?

## Approach

Seven lumped nodes: six panels plus one internal node. Body axes follow LVLH
(+X ram, +Z nadir). Per panel:

```
Q_abs = alpha A [ S max(0, n.s) * sunlit  +  a S F max(0, cos psi) ]  +  eps A q_IR F
Q_out = eps sigma A T^4
```

- `F` is the plate-to-Earth view factor for a plate whose normal is `lam`
  from nadir, `H = r/R`, `sin phi = 1/H`. Fully visible disk: `F = cos(lam)/H^2`.
  Partial visibility uses the closed form from Howell's catalog (case B-22).
- `psi` is the sun angle at the sub-satellite point (simple cosine albedo model).
- Eclipse uses a cylindrical shadow; the analytic eclipse fraction is
  `acos( sqrt(h^2 + 2Rh) / (r cos beta) ) / pi`.
- Adjacent panels are linked by 0.25 W/K joints, each panel to the stack by
  0.12 W/K standoffs plus cavity radiation `sigma eps_int A (T_j^4 - T_i^4)`.

`C dT/dt = Q(t) + G T + radiation(T)` is marched with backward Euler and a
Newton solve per step (analytic Jacobian of the T^4 terms). Whole orbits are
repeated until the start-of-orbit state changes by less than 1 mK.

## Verification (`test_thermal_model.py`)

| Check | Reference |
|---|---|
| Earth view factor, 21 cases at 400-2000 km, 0-120 deg | Monte Carlo cosine-weighted ray casting, within 2e-3 |
| Vertical plate view factor | `(phi - sin phi cos phi)/pi`, 1e-10 |
| Eclipse fraction from the shadow test | closed form above, 1e-4; hand value 0.3777 at beta = 0 |
| Sub-solar point face loads | hand calculation for zenith (full sun) and nadir (albedo + IR) |
| Single-node radiative cooling | `1/T^3 = 1/T0^3 + 3kt/(mc)`, error halves with dt (first order) |
| Isothermal sphere in sunlight | `(alpha S / 4 eps sigma)^(1/4)` = 278.3 K |
| Periodic orbit | energy in = energy out within 0.5%; +Y/-Y identical at beta = 0 |
| High-conductance spacecraft | mean T within 1.5 K of the orbit-averaged single-node balance |
| Time step 360 vs 1440 per orbit | extremes agree within 1 K |

## Run

```
pip install numpy pytest
python3 thermal_model.py
python3 -m pytest -q
```

## Sample output (beta = 0 and 75 deg)

```
--- 500 km circular, beta = 0 deg ---
orbit period 5668.1 s, eclipse fraction 0.379, periodic after 13 orbits
     node    min C   mean C    max C  swing K
       +X    -31.0    -19.6     -6.0     25.0
       -X    -31.0    -19.6     -4.1     26.9
       +Y    -32.0    -21.4     -7.5     24.5
       -Y    -32.0    -21.4     -7.5     24.5
       +Z    -29.1    -19.5     -7.4     21.6
       -Z    -33.7    -16.1     10.4     44.2
 internal    -22.1    -15.7     -9.2     12.9

--- 500 km circular, beta = 75 deg ---
orbit period 5668.1 s, eclipse fraction 0.000, periodic after 8 orbits
     node    min C   mean C    max C  swing K
       -Y     35.5     37.4     39.9      4.4
 internal     22.7     23.9     25.2      2.5
```

The takeaway: at low beta the internal stack sits around -16 C on average and
dips to -22 C, below the 0 C charge limit of a typical Li-ion pack, so the
battery needs a heater or more isolation from the panels. At high beta the
sun-facing -Y panel runs near 40 C with almost no swing, and the stack warms
to about 24 C. The cold case, not the hot case, drives this design.

## Limitations

Lumped panels (no in-plane gradients), a cosine albedo model, mean Earth IR
and albedo instead of hot/cold design cases, and a cylindrical (not conical)
shadow, so penumbra is ignored.
