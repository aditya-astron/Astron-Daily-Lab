# 061 - Sector TID Analysis and Spot-Shield Sizing for an Avionics Box

**Domain:** Radiation shielding

## Problem

A 15-year GEO bus carries a 30 x 20 x 15 cm avionics box (2 mm Al walls)
bolted to the inside of its +X panel. The box holds five cards, and four
parts on them have to stay inside a 2x radiation design margin (RDM):
predicted total ionising dose x 2 <= part capability.

A single "equivalent shielding" number tells you little. Behind the 2 mm wall
alone the curve gives 99 krad. The thinnest line of sight (panel + wall,
0.94 g/cm^2) gives 31 krad, which every part would fail. The real answer
depends on how much mass each part sees in every direction. The questions:

1. What dose does each part actually see, and which fail?
2. Where does that dose come from?
3. What is the lightest fix: spot-shield cans, a wall doubler, or thicker
   box walls?

## Approach

**Mass model.** The bus is a set of axis-aligned boxes and spheres with
densities: honeycomb panels (0.40 g/cm^2), a titanium propellant tank at 50 %
average fill, a battery pack, a wheel/harness block, the box shell and
FR-4 cards. Hollow bodies are an outer solid plus an inner solid of negative
density, so a ray's areal density is just `sum rho_i * chord_i`. Each chord
comes from the slab method (boxes) or the ray-sphere quadratic. All
materials are treated as aluminium-equivalent by mass thickness.

**Sector integration.** For an isotropic environment, the dose at a point
is the solid-angle average of the solid-sphere dose-depth curve:

```
D(p) = (1/4pi) * integral D_sphere(t(p, u)) dOmega
```

Directions come from a Fibonacci lattice (equal-area cells; midpoint rule in
`z`). Section 3 shows it is 10x to 170x more accurate than random sampling at
the same N (10k to 100k directions) in the full model.

**Dose-depth curve.** `D(t) = sum A_k exp(-t/lambda_k)` with three terms
(trapped electrons, harder electrons, protons/bremsstrahlung). The
parameters are **illustrative**, shaped like a 15-year GEO curve (424 krad
at 1 mm, 2.1 krad at 10 mm). For flight work, replace them with
SHIELDOSE-2 or OMERE output. The sum-of-exponentials form also gives an
exact check: between two infinite plates of areal density w, averaging over
directions gives `sum A_k E2(w/lambda_k)`.

**Sizing.** Each fix is a geometric change that gets fully re-traced, and
bisection finds the thickness that puts the part on `capability / RDM`:

- 5-sided Al can (lid + four walls standing on the card, 2.4 cm footprint,
  8 mm tall)
- the same can plus a matching plate under the card
- an Al doubler on the inside of the +X wall
- all six box walls thickened inward

## Verification

`cargo test --release` runs 13 tests:

- E1 against A&S/DLMF tables to 1e-13. E2 against its defining integral on
  both sides of the series/continued-fraction switch.
- Box and sphere chords from inside, outside, parallel and diagonal.
  Cavity subtraction and shell masses are exact.
- Fibonacci lattice: `<z^2>` has an error of exactly `1/(3N^2)`.
- Sector dose against closed forms:
  - centred sphere shell: exact
  - off-centre shell vs 1D quadrature: 1e-9
  - two infinite plates vs E2: 1e-9
  - cubic shell vs 2D face quadrature: 1e-6
- Fibonacci converges, and beats random sampling, in the bus model.
- Spot-can mass against the hand formula. A closed can does better than
  adding w along every ray, and the open can does worse than the closed one.
- Bisection lands on the requirement, and dose decreases with thickness.

`python3 verify/crosscheck.py` recomputes the references with scipy
(`expn`, `quad`, `dblquad`) and re-traces the bus model with independent
vectorised numpy code. All values agree with the Rust output to its printed
precision.

## Run

```
cargo run --release              # design study (~7 s)
cargo test --release             # verification
python3 verify/crosscheck.py     # scipy/numpy cross-check (needs numpy, scipy)
```

## Sample output (abridged)

```
2. Sector integrator against closed forms (N = 200000 Fibonacci directions)
   sphere shell, centre        43.363633  exact  43.363633  rel +2.2e-13
   sphere shell, d = 0.6 R     32.020346  quad   32.020346  rel -4.3e-12
   two infinite plates         31.621533  E2     31.621533  rel -2.7e-11
   cubic shell, centre         68.789136  quad   68.789132  rel +5.8e-8

3. Convergence at U2 SRAM (card 5, centre) in the full model
        N    Fibonacci err    random err   random 1-sigma
    10000       -6.56e-4      -7.06e-3        1.14e-2
   100000       +1.83e-5      +3.07e-3        3.62e-3

4. Baseline sector doses, requirement D * 2 <= capability
   part                       D krad  eq.sphere  t_min  D(t_min)  t<1 g/cm^2   cap  margin
   U1 FPGA (card 5, corner)    6.89   6.15 mm  0.940      31.4        5.2 %   100  14.51 pass
   U2 SRAM (card 5, centre)    7.31   6.02 mm  0.940      31.4        5.1 %    10   1.37 FAIL
   U3 ADC  (card 3, centre)    5.11   6.81 mm  0.940      31.4        3.4 %    20   3.91 pass
   U4 LDO  (card 1, corner)    5.46   6.65 mm  0.940      31.4        4.5 %    10   1.83 FAIL

5. Where the dose comes from (share of D by avionics-box face crossed)
   U2 SRAM:  +X (panel) 52.2% (21.2% of 4pi) ...
   U4 LDO :  +X (panel)  1.5% ( 2.8% of 4pi)  ...  -Z 36.9% (23.6% of 4pi)

6. Mitigation trade for the failing parts
   5-sided can on U2 SRAM: wall 0.69 mm,    2.7 g
   5-sided can on U4 LDO : wall 0.17 mm,    0.6 g
   6-sided can + back plate on U2 SRAM: wall 0.59 mm,    3.3 g
   +X doubler for U4 LDO : not enough even at 30 mm (5.38 krad, needs 5.00)
   +X doubler 1.39 mm: 107.5 g (does not close the design)
   uniform wall 2.00 -> 2.63 mm: 439.9 g

7. Chosen design: 5-sided cans, all parts re-evaluated together
   U2 SRAM (card 5, centre)    7.31    5.00   2.00 pass
   U4 LDO  (card 1, corner)    5.46    5.00   2.00 pass
   added mass 3.3 g vs 439.9 g for uniform walls (134x)
```

## Takeaways

- Sector analysis gives 5-7 krad. The thinnest-path bound gives 31 krad.
  The box looks like 6-7 mm of equivalent Al even though its wall is 2 mm,
  because the cards, the other components and the tank fill most of the sky.
- Dose follows the thinnest solid angle, not the nearest wall. Half of the
  SRAM's dose enters through the +X panel, which covers only 21 % of its sky.
  The LDO sits behind four cards from that side, so it gets 1.5 % from +X and
  most of its dose from the -Z and -Y walls.
- A +X doubler (107 g) fixes the SRAM and leaves the LDO failing even at
  30 mm. Uniform walls close the design at 440 g. Two spot cans close it at
  3.3 g, so local shielding is about two orders of magnitude lighter.
- A backside plate adds mass and saves little. The card and the
  rest of the stack already shield the back hemisphere, so the open can is
  lighter. The computed 0.17 mm LDO can would in practice be set by minimum
  gauge (about 0.5 mm), which adds margin.
