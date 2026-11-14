# 037 - Water Recovery Sorbent Bed: Breakthrough and LUB Scale-Up

**Domain:** ECLSS (water recovery)

## Problem

The polishing stage of a 4-crew water processor passes 600 mL/h (14.4 kg/day)
of distillate through a 10 cm ID bed of granular activated carbon. The feed
carries a model organic contaminant at 20 mg/L with a Langmuir isotherm
(`qm` = 100 mg/g, `b` = 0.2 L/mg, so `b c0` = 4 and separation factor
`R = 1/(1 + b c0)` = 0.2). Uptake is intraparticle-limited, `k` = 0.005 1/h.
The product limit is 1 % of the feed. How long must the bed be to last a year,
and how much does the mass-transfer zone cost compared with an equilibrium
estimate?

## Approach

**Model.** Plug flow with a linear-driving-force uptake:

    eps dc/dt + eps u dc/dz + rho_b dq/dt = 0
    dq/dt = k (q*(c) - q),      q* = qm b c / (1 + b c)

**Characteristic grid.** With the retarded time `tau = t - z/u` the PDE splits
into two ODEs, one along each characteristic:

    dc/dz  |tau = -(rho_b k / (eps u)) (q* - q)
    dq/dtau|z   =  k (q* - q)

with `c(0, tau) = c0` and `q(z, 0) = 0`. Both are integrated with the
trapezoidal rule, so each grid node is a single monotone equation in `c`
(Newton, `g' > 0`). No numerical dispersion is introduced, the scheme is
second order, and mass is conserved to round-off. The trapezoid is A- but not
L-stable, so the solver refuses a z step holding more than ~2 clean-bed
transfer units per cell.

**Constant pattern.** In a developed wave `q/q0 = c/c0 = x`, and

    dx/dt = k (1 - R) x (1 - x) / (R + (1 - R) x)

integrates in closed form. Fixing the constant with the mass balance
`int (x - H(t - t_st)) dt = 0` gives

    tau(x) = tau_st - 1/k + [R ln x - ln(1 - x)] / (k (1 - R))
    tau_st = rho_b q0 L / (eps u c0)

**LUB scale-up.** The length of unused bed `LUB = L (1 - t_b / t_st)` from a
pilot run is added to the equilibrium length for the target life,
`L = w t_life + LUB`, `w = L_pilot / t_st` the wave speed. The scaled bed is
then simulated directly.

## Verification (`make test`)

- Anzelius/Thomas function `J(xi, eta) = 1 - int_0^xi e^(-eta-s) I0(2 sqrt(eta s)) ds`:
  `J(xi, 0) = e^-xi` and the identity `J(x,y) + J(y,x) = 1 + e^(-x-y) I0(2 sqrt(xy))`
  to 4e-16.
- Linear isotherm (xi = 20 transfer units): solver against the exact Thomas
  solution, max error 4.0e-6 on 400x1600, observed order 1.999.
- First moment of the breakthrough curve equals `tau_st` (linear: 6e-6 h of 40 h;
  Langmuir: 2e-10 h of 9000 h).
- Langmuir bed 15 MTZ lengths long: breakthrough times at 1, 10, 50, 90, 99 %
  match the constant-pattern formula within 0.32 h on a 138 h wide zone.
- A z grid with 11 transfer units per cell is rejected.

Hand check of the design inputs: `q0 = 100 * 4/5 = 80 mg/g`,
`u = 600 / (78.54 * 0.4) = 19.10 cm/h`,
`t_st = 0.45 * 80 * 20 / (0.4 * 19.10 * 0.02) = 4712 h`, and the
constant-pattern LUB `(1 + 1.139)/k * 20 / 4712 = 1.816 cm`.

## How to run

    make test    # verification
    make run     # pilot bed, LUB scale-up, simulated check

Needs a C++17 compiler (uses `std::cyl_bessel_i`).

## Sample output

    Pilot bed L = 20.0 cm (0.71 kg carbon)
      stoichiometric time     4713.4 h  (196.4 days)
      1% breakthrough         4285.7 h  (178.6 days)
      50% / 99%               4652.1 h / 5664.2 h
      utilisation at limit     0.909
      LUB                      1.815 cm
      constant-pattern LUB     1.815 cm

    Sizing for 365 days at c/c0 <= 0.01
      equilibrium length (no MTZ)    37.18 cm
      LUB scale-up length            38.99 cm
      simulated life of that bed     365.1 days (target 365, error +0.02 %)
      carbon mass                     1.38 kg
      equilibrium-sized bed lasts    347.3 days (17.7 days short)

The 20 cm pilot already holds a fully developed constant-pattern front (its
LUB equals the asymptote), which is why the LUB shortcut lands within 0.02 %.
Sizing on equilibrium capacity alone would leave the crew 18 days short of a
year. The long tail (99 % takes another 40 days after 50 %) is the favourable
isotherm's signature: `R ln x` keeps the leading edge sharp, `-ln(1 - x)`
drags the trailing edge out.
