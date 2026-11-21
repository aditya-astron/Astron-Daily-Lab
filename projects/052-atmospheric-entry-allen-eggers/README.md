# 052 - Ballistic and lifting Earth entry: Allen-Eggers, equilibrium glide and Sutton-Graves heating

A capsule or probe coming back from LEO (7.8 km/s) or the Moon (11 km/s)
has to shed its kinetic energy within a peak-g limit for crew or hardware,
and within a peak heat flux and total heat load its thermal protection can
take. Those three numbers pull in different directions with entry angle and
ballistic coefficient. This project integrates planar point-mass entry,
checks the integrator against the classical closed forms to solver
precision, then measures where those closed forms stop being trustworthy and
runs the trade.

## Approach

**Dynamics.** Non-rotating spherical Earth, isothermal exponential
atmosphere (rho0 = 1.225 kg/m^3, H = 7.2 km), inverse-square gravity:

```
dV/dt     = -D/m - g sin(gamma)            D/m = rho V^2 / (2 beta),  beta = m/(CD A)
dgamma/dt = L/(m V) - (g/V - V/r) cos(gamma)
dh/dt     = V sin(gamma)
ds/dt     = V cos(gamma) Re/r
dQ/dt     = k sqrt(rho/rn) V^3             Sutton-Graves, k = 1.7415e-4
```

DOP853 at rtol 1e-11 with dense output; peaks are refined on the
interpolant and the heat load is carried as a state, so neither depends on
output sampling. Events stop the run at the ground, at a cutoff speed, or on
a skip-out back through the interface.

**Allen-Eggers (gravity and curvature off, gamma frozen).** With s = sin|gamma_e|:

```
V(rho)     = Ve exp(-rho H / (2 beta s))
peak decel = Ve^2 s / (2 e H)              at rho = beta s / H,      V = Ve / sqrt(e)
peak q     = k Ve^3 sqrt(beta s / (3 e H rn)) at rho = beta s / (3 H), V = Ve e^(-1/6)
heat load  = k Ve^2 sqrt(pi beta H / (rn s)) [erf(sqrt(rho_g H/(beta s))) - erf(sqrt(rho_e H/(beta s)))]
```

The heat-load integrand goes as rho^(-1/2), so the thin air above the 120 km
interface is not negligible: ignoring it is a 0.5 % error at 5 deg. Both the
rho_e term and the matching rho = 0 speed are carried explicitly.

**Equilibrium glide.** Lift balances weight less centrifugal relief,
L/m = g (1 - V^2/Vc^2), which gives

```
decel = g (1 - V^2/Vc^2) / (L/D)
range = (r L/D / 2) ln[(1 - Vf^2/Vc^2) / (1 - V0^2/Vc^2)]
sin(gamma_eq) = -2 H g / (V^2 L/D)
```

The glide altitude is a fixed-point solve because g and r move with h.

## Results

- Integrator vs. Allen-Eggers with gravity and curvature off: peak g, peak
  flux and heat load agree to about 1e-11, V(rho) to 1e-9 along the path.
- With the full model, steep entries run about 3 % above the closed-form
  g and 4 % above the flux: the kinetic energy gravity adds on the way down.
  Below about 5 deg the closed forms fail badly (peak g off by 2x at 1.5 deg)
  because gravity and curvature steepen the path into the dense layer.
- At 7.8 km/s peak g is independent of beta; flux and load scale as
  sqrt(beta sin gamma) and sqrt(beta / sin gamma). Steep is a short, intense
  pulse; shallow is a long soak.
- At 11 km/s a 2 deg ballistic entry skips back out, and at 6 deg high beta
  trades lower peak g for much larger flux and load.
- The equilibrium-glide range formula is within 1-2 % at L/D 1-2 and 8 %
  out for a capsule at L/D 0.3.

## Run

```
pip install numpy scipy pytest
python run_entry.py       # full study, about 1 s
python -m pytest -q       # 12 checks
```

## Sample output (excerpt; full run in sample_output.txt)

```
2. Full model (gravity + curvature) vs. Allen-Eggers, Ve = 7.8 km/s, beta = 300, rn = 0.5
    gamma  g full    g AE     err  q full    q AE     err h_q full  h_q AE
     -1.5    8.52    4.15 +105.4%   1.614   1.352  +19.4%     54.7    58.5
     -5.0   15.55   13.81  +12.6%   2.592   2.466   +5.1%     49.3    49.8
    -20.0   56.18   54.21   +3.6%   5.083   4.886   +4.0%     39.9    40.0
    -90.0  163.76  158.49   +3.3%   8.714   8.354   +4.3%     32.2    32.3

4. Equilibrium glide from 7.5 km/s to 1.5 km/s, beta = 400 kg/m^2
    L/D  h0 km gamma0 deg  range km closed km     err  peak g   eq g
    0.3   74.2     -0.469      2441      2261   +8.0%    2.46   3.28
    1.0   83.0     -0.140      7727      7580   +1.9%    1.24   1.33
    2.0   88.0     -0.070     15402     15211   +1.3%    1.04   1.05
```

## Tests

`test_entry.py` checks the integrator against every Allen-Eggers closed form
at two speeds and three angles, the erf heat-load expression against direct
quadrature, hand-computed peak g and altitude for a 20 deg entry, the steep-
entry gravity bias against 2 g dh / Ve^2, beta-independence of peak g at
circular speed, skip-out at 11 km/s and 2 deg, and the equilibrium glide
state, altitude tracking and range.

## Limits

Single-scale-height atmosphere, no Earth rotation, constant aero
coefficients, convective stagnation heating only (no radiative term, which
matters above about 10 km/s for large noses).
