# 014 - Oxygen loop closure: PEM electrolysis and Sabatier CO2 reduction

## Problem

A 4-crew habitat makes its oxygen by electrolysing water and reduces the
exhaled CO2 with the by-product hydrogen in a Sabatier reactor. How much of
the electrolyser's water comes back, how much makeup water has to be
launched per year, and what do the stack and reactor cost in power and heat?

## Approach

**Thermochemistry.** NIST Shomate fits for H2, O2, CO2, CH4, H2O(g) and
H2O(l) give H, S and G at any temperature, so the equilibrium constant is
`ln K = -dG(T) / RT` with no constant-dH shortcut.

**Sabatier equilibrium.** `CO2 + 4 H2 = CH4 + 2 H2O` at 1 bar:

    Kp = y_CH4 y_H2O^2 / (y_CO2 y_H2^4) * (P/P0)^-2

solved for the extent `xi` with Brent's method (ln Q is monotonic in xi,
so the root is unique). RWGS to CO is neglected.

**Loop balance.** Electrolysis gives `n_H2 = 2 n_O2`; the crew returns
`n_CO2 = RQ n_O2` with RQ ~ 0.9. The feed H2/CO2 is ~2.2, far below 4, so
hydrogen limits and

    closure = 2 xi / (2 n_O2) <= 1/2

Half the hydrogen leaves as methane: the loop can never recover more than
50% of the water, whatever the reactor does.

**Electrolyser.** Per cell, `V = E_rev(T) + (RT/aF) asinh(i/2i0) + i ASR`;
stack current from Faraday's law (4 e- per O2), waste heat
`(V - V_tn) I`.

## Verification

`test_o2_loop.py` checks:

- Shomate Cp and S against JANAF at 298 K and 1000 K for every species
- Sabatier dH298 = -165.0 kJ, dG298 = -113.5 kJ; E_rev = 1.229 V,
  V_tn = 1.481 V at 25 C, E_rev = 1.184 V at 80 C
- Gibbs-Helmholtz: numeric d lnK/dT equals dH/RT^2
- ln K at 300 C against an independent constant-dCp Kirchhoff hand calc
- solved extent reproduces K exactly and conserves C, H, O
- conversion falls with T, rises with P; 4:1 feed conversion inside the
  published 1 atm ranges at 300 C and 500 C
- daily water (3.783 kg) and stack current (39.09 A) by hand, closure < 50%
  everywhere and -> 50% when cold, stack power = V_tn I + waste heat

## Run

    pip install scipy
    python3 test_o2_loop.py
    python3 o2_loop.py

## Sample output

```
Crew 4: O2 3.36 kg/d, CO2 4.16 kg/d, RQ 0.900
PEM cell at 80 C: E_rev 1.1831 V, V_tn 1.4722 V

Cell polarisation
  i A/cm2   V cell   eff(HHV)
     0.10   1.6216     0.908
     0.50   1.7426     0.845
     1.00   1.8536     0.794
     1.50   1.9560     0.753
     2.00   2.0547     0.716

Stoichiometric 4:1 feed, 1 bar: equilibrium CO2 conversion
   200 C  ln K   20.64  X_CO2 0.9891
   300 C  ln K   12.85  X_CO2 0.9492
   400 C  ln K    7.24  X_CO2 0.8504
   500 C  ln K    3.00  X_CO2 0.6819
   600 C  ln K   -0.32  X_CO2 0.4784
   700 C  ln K   -3.00  X_CO2 0.3015

Closed loop, crew-limited H2 feed, 1 bar
  T_sab C  H2 conv  closure  makeup kg/d  CH4 kg/d  Sab heat W
      250   0.9880   0.4940        1.914     0.832       105.6
      300   0.9713   0.4857        1.946     0.818       105.0
      350   0.9406   0.4703        2.004     0.792       102.8
      400   0.8911   0.4455        2.098     0.751        98.4
      450   0.8203   0.4102        2.232     0.691        91.5
      500   0.7300   0.3650        2.403     0.615        82.1
      600   0.5192   0.2596        2.801     0.437        59.3

Design point: Sabatier 300 C, PEM stack 12 cells x 40 cm2 at 80 C
  H2/CO2 feed ratio      2.222  (4.0 needed for full CO2 reduction)
  water electrolysed     3.783 kg/d
  water recovered        1.837 kg/d  (48.57%, ceiling 50%)
  makeup water           1.946 kg/d  = 710 kg/yr
  CO2 vented             1.916 kg/d
  stack current          39.09 A  (0.977 A/cm2, 1.849 V/cell)
  stack power            867.2 W  (6.19 kWh/kg O2)
  stack waste heat       176.7 W
  Sabatier heat release  105.0 W
```

The reactor runs close to the ceiling at 250-350 C; the remaining ~1.9 kg/d
(about 700 kg/yr for four crew) is the price of venting methane. Pyrolysing
the CH4 to recover its hydrogen, or a Bosch reactor, is what it takes to
break the 50% limit.
