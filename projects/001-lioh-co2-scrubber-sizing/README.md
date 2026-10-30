# 001 - LiOH CO2 Scrubber Sizing for a 4-Crew 30-Day Mission

**Domain:** ECLSS

## Problem

A sealed 20 m^3 cabin holds 4 crew for 30 days. CO2 is removed by expendable
lithium hydroxide canisters. How many canisters must fly, how often are they
swapped, and what side loads (water, heat) does the scrubber put on the cabin?

A pure stoichiometric budget overstates usable capacity, because a bed stops
removing CO2 well before it is fully converted. Here the swap point is set by
the cabin ppCO2 limit, and the usable capacity falls out of the simulation.

## Approach

Chemistry: `2 LiOH + CO2 -> Li2CO3 + H2O`

- LiOH per kg CO2: `2 M_LiOH / M_CO2 = 1.088 kg/kg`
- Product water: `M_H2O / M_CO2 = 0.409 kg/kg`
- Reaction enthalpy from formation enthalpies: -94.5 kJ/mol CO2 (water as vapour)

Well-mixed cabin with one active canister (state: cabin CO2 concentration `C`,
CO2 held in canister `A`):

```
dC/dt = (G - eta(u) Q C) / V
dA/dt = eta(u) Q C,          u = A / A_max
eta(u) = eta0 / (1 + exp((u - u50)/w))   (normalised so eta(0) = eta0)
```

`G` is crew CO2 generation (1.04 kg/crew-day, NASA BVAD nominal), `Q` the fan
flow through the bed. A fresh bed holds the cabin at `C = G / (eta0 Q)`. The
canister is swapped when ppCO2 hits 3.0 mmHg. Integrated with LSODA and a
terminal event at the limit.

Breakthrough parameters (`u50`, `w`) are representative, not from a specific
flight canister; change them in `Canister` to match vendor data.

## Run

```
pip install numpy scipy pytest
python3 lioh_sizing.py
python3 -m pytest -q
```

## Sample output

```
Mission: 4 crew, 30 days, cabin 20.0 m^3
CO2 generated                :   124.80 kg (4.16 kg/day)
LiOH per kg CO2 (theory)     :   1.0882 kg/kg
LiOH, 100% utilisation       :   135.81 kg
Fresh-bed steady ppCO2       :    2.055 mmHg
Swap limit                   :    3.000 mmHg
Canister CO2 capacity        :    3.676 kg (4.0 kg LiOH)
Mean utilisation at swap     :    0.880
Mean canister life           :    18.68 h
Canisters consumed (sim)     :       39
Canisters to fly (+10% spares):       43
LiOH flown                   :    172.0 kg
Product water released       :    51.06 kg
Mean reaction heat load      :    103.3 W
Peak / mean ppCO2            : 3.000 / 2.220 mmHg
Mass balance error           : 3.89e-13 kg
```

Hand check: canister life = 0.88 x 3.676 kg / 0.1733 kg/h = 18.7 h, and
124.8 kg / (0.88 x 3.676 kg) = 38.6, so 39 canisters. The tests cover these
checks plus stoichiometry, enthalpy, steady state and mass conservation.
