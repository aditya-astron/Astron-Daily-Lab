# 045 - Multilayer insulation heat leak, layer by layer

## Problem

How much heat leaks through a multilayer insulation (MLI) blanket, what does
each shield run at, and how many shields should go into a blanket of fixed
thickness? Lumped thermal models carry MLI as one effective emittance eps*;
this project builds that number from the stack itself (radiation between
shields, conduction through the spacer net, residual gas), solves it for
fixed and space-facing boundaries, and checks it against exact results and
the Lockheed empirical equation.

## Approach

Node 0 is the inner layer held at the wall temperature, nodes 1..n are
shields, node n+1 is the cold boundary. Each gap carries

    q = sigma (Ta^4 - Tb^4) / (1/ea + 1/eb - 1)     grey parallel plates
      + (hs/2)(Ta^2 - Tb^2),   hs = Cs Nbar^2.63    spacer net
      + Cg P (Ta^0.52 - Tb^0.52)                    free-molecular gas

The spacer term is a conductance hs * T_mean per gap. Summed over n+1 gaps it
telescopes to exactly the Lockheed solid term Cs Nbar^2.63 (Th^2 - Tc^2) /
(2(n+1)), so the layer model and the empirical fit share one conduction law
(Cs = 8.95e-8, Cr = 5.39e-10, Cg = 1.46e4; Nbar in layers/cm, P in torr).

**Solvers.**

- Newton on the nodal balances q_(i-1) - q_i = 0 with an analytic tridiagonal
  Jacobian (Thomas sweep). A radiating outer cover adds one equation,
  q_n + alpha S = eps_o sigma T^4.
- Shooting, as an independent check: guess q, march gap by gap solving for
  the next temperature, bisect q on the far boundary condition.
- Closed form for uniform emissivity: every gap flux is F(Ta) - F(Tb) for a
  single potential F(T) = sigma T^4/(2/eps - 1) + hs T^2/2 + Cg P T^0.52, so
  q = (F(Th) - F(Tc))/(n+1) and T_k = F^-1(F(Th) - k q).

**Shield count at fixed thickness.** With thickness t, Nbar = n/t. Radiation
falls as 1/(n+1) while conduction grows as n^2.63/(n+1), so q(n) has an
interior minimum, found by golden section on the closed form and confirmed by
integer brute force and by the Newton layer model.

## Checks (all in `mli_test.go`)

| Check | Result |
|---|---|
| Pure radiation vs sigma dT^4 / ((n+1)(2/eps-1)), n = 1..60 | rel err < 1e-12, T^4 linear |
| Radiation off vs Lockheed solid term | exact to 1e-12 |
| Newton vs closed form vs shooting (uniform stack) | 1e-13 K |
| Non-uniform stack (paint, crinkled shields): Newton vs shooting, equal gap fluxes | 1e-14 rel |
| Radiating cover: cover energy balance, Newton vs shooting | 1e-10 W/m^2 |
| Zero leak when absorbed sun puts the cover at wall temperature | S = 0.75 sigma 293^4 / 0.40 = 783.6 W/m^2 |
| Fixed-thickness optimum: integer argmin adjacent to continuous, df/dn = 0 | n* = 14.7, best n = 15 |
| Ideal-shield model vs Lockheed fit at 300/77 K, 20 layers | within 4.3 % |

Hand check for case 1: sigma/(2/0.03 - 1) = 8.635e-10, 293^4 - 150^4 = 6.864e9,
over 21 gaps gives 0.2822 W/m^2.

## Run

    go run .
    go test -v .

## Sample output

    == 2. Lab boundary 300 K -> 77 K, 20 shields, Nbar 20 /cm, 1e-6 torr ==
    Newton      q = 0.822767 W/m^2 (5 iters)
    closed form q = 0.822767  (F(Th)-F(Tc))/(n+1), max |dT| 5.68e-14 K
    shooting    q = 0.822767  max |dT| 4.26e-13 K
    model    : rad 0.3429  solid 0.4731  gas 0.0068  total 0.8228 W/m^2  eps* 0.00180
    Lockheed : rad 0.3085  solid 0.4731  gas 0.0072  total 0.7888 W/m^2  eps* 0.00172

    == 4. Flight blanket on a 293 K wall: 18 shields, cover alpha 0.40 eps 0.75 ==
    eclipse     cover   60.6 K  q into blanket  +0.5750 W/m^2  (shoot  +0.5750)  eps* vs cover 0.00138
    sun 45 deg  cover  308.4 K  q into blanket  -0.1034 W/m^2  (shoot  -0.1034)  eps* vs cover 0.00109
    sun normal  cover  336.3 K  q into blanket  -0.3261 W/m^2  (shoot  -0.3261)  eps* vs cover 0.00106
    zero-leak absorbed-sun flux 783.58 W/m^2 (cover at 293.00 K, radiative check 293.00 K)

    == 5. Fixed 1.0 cm blanket, 293 K -> 150 K: how many shields? ==
    continuous optimum n* = 14.68 (q 0.5949), integer brute force n = 15 (q 0.5951 W/m^2)

## Notes

- Conduction through the net already exceeds radiation at 20 layers/cm; past
  about 15 shields per cm, adding layers makes the blanket worse.
- The radiation share sits on the cold side: the last shields run far apart
  in temperature (109 K to 77 K) because T^4 differences are small there.
- These are ideal-blanket numbers (eps* ~ 1e-3). Flight blankets with seams,
  vents, pins and edge effects typically show eps* of 0.01 to 0.03, so use
  this model for the stack and add a measured degradation for the installation.
- The eclipse cover temperature is the steady-state limit; a real cover is
  still cooling when the next sunrise arrives. Cover optics are assumed values.
