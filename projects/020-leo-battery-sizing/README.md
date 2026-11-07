# 020 - LEO Li-ion Battery Sizing with a 2RC Equivalent-Circuit Model

**Domain:** Spacecraft power systems

## Problem

A small LEO spacecraft at 500 km draws a constant 180 W from a regulated bus.
Its battery is built from 3.2 Ah 18650-class NMC cells, 8 in series. How many
parallel strings are needed to keep depth of discharge at or below 25% at the
worst-case beta angle, what does the pack voltage do over an orbit, how much
of the sunlit arc is spent tapering on the charge-voltage limit, and how much
array power does it take before the battery is limited by its own charge
acceptance rather than by power?

The usual hand sizing (eclipse energy / nominal voltage / DoD) answers the
first question approximately. This project checks it against a cell model
that has real resistance, relaxation and a CC-CV charger, run to periodic
steady state.

## Approach

**Eclipse geometry.** Cylindrical shadow, circular orbit, Sun at beta above
the orbit plane:

```
cos(pi f) = sqrt(h^2 + 2 R h) / (r cos beta),   no eclipse once sin(beta) >= R/r
```

At 500 km and beta = 0 that is 35.75 min of a 94.62 min orbit. Eclipses stop
above beta = 68.02 deg.

**Cell model.** Thevenin circuit with two RC branches (I > 0 is discharge):

```
V_t = OCV(SOC) - V1 - V2 - I R0
dVk/dt = -Vk/(Rk Ck) + I/Ck          tau1 = 20 s, tau2 = 600 s
dSOC/dt = -I/Q
```

The regulator fixes current at the start of each step, so the RC states are
advanced with the exact zero-order-hold solution
`Vk+ = a Vk + Rk (1 - a) I`, `a = exp(-dt/tau)`.

**Power flow.**
- Eclipse: the battery supplies `P_load / eta_BDR` at constant power. The
  current is the stable root of `P = (E - I R0) I`, `E = OCV - V1 - V2`; no
  root means voltage collapse.
- Sunlight: the surplus `(P_array - P_load) eta_BCR` is offered to the
  battery. Charge current is the smallest of the power-limited current, the
  C/2 current limit and `(V_max - E)/R0`, which keeps the terminal at
  4.10 V/cell. That last term produces the CV taper naturally. Whatever the
  battery does not take is shunted.

**Steady state.** Orbits are repeated from a full, rested pack until the
eclipse-entry state repeats to 1e-10.

**Energy accounting.** Terminal energy in minus out must equal the change in
chemical energy `Q * integral OCV dSOC` (exact for the piecewise-linear OCV),
plus the change in RC capacitor energy, plus Joule loss in R0 and the RC
resistors (integrated exactly over each step). Terminal power is integrated
the way a power monitor would, sampled and held, so the residual should fall
as O(dt). It does.

**Array threshold.** Bisection on array power for the smallest value whose
steady-state eclipse-entry SOC is within 0.2% of the value reached with
unlimited power. Compared against the plain energy balance
`P_sa = P_load + E_batt / (eta_rt eta_BCR T_day)`.

## Verification (`tests/verify.rs`)

- Orbit period and eclipse fraction against hand values; eclipse vanishes at
  the critical beta and decreases monotonically with beta.
- RC branches match `I R (1 - exp(-t/tau))` to 1e-12 for step sizes from
  0.1 s to 45 s, and relax as `exp(-t/tau)` at zero current; Coulomb counting
  is exact.
- Closed-form OCV integral against midpoint quadrature.
- Constant-power current roots deliver exactly the requested power, and none
  exists beyond `E^2/(4 R0)`.
- Lossless cell: end-of-eclipse SOC matches `Q (W(s0) - W(s1)) = P t_e`
  solved independently.
- Steady state is periodic, never exceeds 4.10 V/cell or C/2.
- Energy balance residual below 1e-4 and halves when dt halves.
- Sizing returns 5 strings, 4 strings fails DoD, and the hand estimate agrees.
- Array threshold sits above the energy-balance estimate and 3% less power
  fails.

## Run

```
cargo run --release
cargo test --release
```

## Sample output

```
=== Sizing (beta = 0, load 180 W, DoD <= 25%) ===
eclipse battery energy   112.91 Wh
hand capacity             15.26 Ah  -> 5 strings of 3.2 Ah
  np  hand DoD[%]  DoD[%]  Vmin[V]  Vmax[V]  Ichg[A]  t_CV[min]  round-trip
   1  battery runs down: eclipse draw exceeds capacity or C/2 recharge
   2  battery runs down: eclipse draw exceeds capacity or C/2 recharge
   3        39.73   38.85    28.75    32.80     4.80       22.2     0.9497
   4        29.80   28.45    29.92    32.80     6.40       36.0     0.9603
   5        23.84   22.49    30.55    32.80     6.71       36.8     0.9685
   6        19.87   18.60    30.93    32.80     6.65       36.9     0.9739
   7        17.03   15.86    31.19    32.80     6.61       37.0     0.9778
simulated sizing       5 strings (8s5p, 16.0 Ah, 474 Wh nameplate)

  energy-balance residual vs step size (first order expected)
    dt =  4.0 s   residual 1.155e-4
    dt =  2.0 s   residual 5.763e-5
    dt =  1.0 s   residual 2.878e-5
    dt =  0.5 s   residual 1.438e-5

=== Array power needed to recharge every orbit ===
acceptance-limited entry SOC         0.9131
bisected (entry SOC within 0.2%)      393.2 W
energy-balance estimate               305.1 W
ratio                                 1.289
```

## Takeaways

- The hand sizing picks the right string count, and the simulated DoD comes
  in about 1.4 points lower. The pack spends the eclipse in the upper quarter
  of its curve, averaging about 3.9 V/cell at the terminals rather than the
  3.70 V nominal.
- 2p fails even though its capacity would cover the eclipse. A C/2 charge
  limit cannot put 3.9 Ah back in a 59 min sunlit arc, so the pack runs down
  over several orbits. Charge acceptance constrains LEO packs as much as
  capacity does.
- The pure energy balance underestimates the array by about 30%. At 400 W the
  pack spends 37 of the 59 sunlit minutes on the voltage limit, where it
  cannot use extra power, so the bulk charge has to happen faster.
- Even with unlimited array power the pack enters eclipse at 0.913 SOC, short
  of the 0.920 open-circuit full point. The 600 s diffusion branch keeps the
  terminal voltage high, so the CV taper ends before the cell is truly full.
