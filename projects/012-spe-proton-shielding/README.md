# 012 - Solar proton event dose behind aluminium vs polyethylene shielding

## Problem

A crew shelter has to ride out a large solar particle event. Using the
August 1972 event as the design case, how much dose gets through a slab wall,
and how much mass does hydrogen-rich polyethylene (HDPE) save over aluminium
for the same dose limit?

## Approach

**Stopping power.** Bethe formula for protons, with the exact `Tmax` and
ICRU mean excitation energies (water 75 eV, Al 166 eV, HDPE 57.4 eV):

    S = K (Z/A) / b^2 [ 0.5 ln(2 me c^2 b^2 g^2 Tmax / I^2) - b^2 ]

Per gram, HDPE stops protons harder than Al because it has more electrons
per gram (Z/A 0.570 vs 0.482) and a lower `I`.

**Range.** CSDA range `R(E) = integral dE / S` tabulated on a log grid from
1 to 3000 MeV, inverted by log-log interpolation.

**Transport.** Straight-ahead, continuous slowing down. A proton leaving with
`E_out` entered with `E_in = E(R(E_out) + x)`; conservation of protons gives
`phi_out = phi_in(E_in) S(E_in) / S(E_out)`, chained across layers.

**Dose.** Thin tissue (water) layer just behind the wall,
`D = 1.602e-10 * integral phi_out S_water dE` Gy. For an isotropic event the
point sees `t / mu` of wall from each direction of the exposed hemisphere:

    D_iso(t) = 1/2 integral_0^1 D_normal(t / mu) dmu      (48-point Gauss-Legendre)

**Spectrum.** King (1974) fit, `Phi(>E) = 7.9e9 exp((30 - E)/26.5)` cm^-2,
cut at 10 MeV.

**Limits.** No nuclear interactions or secondary neutrons, and no shell
corrections (the 10 MeV water range is 1.4% short of PSTAR). Nuclear
secondaries are produced more in Al than in HDPE but also add dose behind
both, so the HDPE advantage beyond ~10 g/cm^2 is optimistic here; the
thin-to-moderate wall numbers are the useful ones.

## Verification

`test_shield.py` checks:

- water stopping power and CSDA range against NIST PSTAR at 10, 100, 200 MeV
- exact range-table inversion and layer additivity (3 + 4 = 7 g/cm^2)
- proton conservation: transmitted count equals `Phi_in(> E(R(1 MeV) + t))`
- unshielded dose against an independent `scipy.integrate.quad`
- `D_iso(0) = D_free / 2`, continuity as `t -> 0`, quadrature convergence
- dose monotonic in thickness, HDPE below Al at every areal density

## Run

    pip install numpy scipy
    python3 test_shield.py
    python3 shield.py

## Sample output

```
King (1974) fit, August 1972 SPE, isotropic, 2-pi slab exposure
Thin-tissue dose behind slab, Gy (no nuclear attenuation)

 g/cm^2     aluminium  polyethylene         water   HDPE/Al
    0.0         29.13         29.13         29.13     1.000
    0.5        14.639        11.018        11.691     0.753
    1.0         8.289         5.768         6.221     0.696
    2.0         3.881         2.425         2.674     0.625
    5.0         0.927         0.467         0.539     0.504
   10.0         0.192         0.076         0.092     0.394
   20.0         0.021         0.006         0.007     0.271
   30.0         0.004         0.001         0.001     0.200

Areal density to hold this event to the dose, g/cm^2:
  2.00 Gy   Al   3.20 ( 1.19 cm)   HDPE   2.28 ( 2.40 cm)   mass saving 29.0%
  1.00 Gy   Al   4.80 ( 1.78 cm)   HDPE   3.43 ( 3.61 cm)   mass saving 28.5%
  0.50 Gy   Al   6.75 ( 2.50 cm)   HDPE   4.85 ( 5.10 cm)   mass saving 28.2%
  0.25 Gy   Al   9.04 ( 3.35 cm)   HDPE   6.52 ( 6.87 cm)   mass saving 27.9%
```

For this event HDPE meets the same dose limit with about 28% less wall mass,
though it needs roughly twice the thickness of Al.
