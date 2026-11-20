# 048 - ADC anti-alias filter: Chebyshev LC ladder vs. Sallen-Key cascade

A 40 MS/s telemetry ADC sits behind a 50 ohm IF chain. Everything from
fs - 8 MHz = 32 MHz upward folds into the 0-8 MHz band, so the filter must give
0.5 dB ripple to 8 MHz and at least 60 dB at 32 MHz. Production acceptance
allows 1.0 dB peak-to-peak ripple after tolerances. The same Chebyshev
response is built two ways, a doubly terminated LC ladder and a unity-gain
Sallen-Key cascade, and the question is which one survives real parts.

## Approach

**Order.** `n = ceil(acosh(sqrt((10^(As/10)-1)/(10^(Ap/10)-1))) / acosh(ws/wp))`
gives n = 5 for Chebyshev (74.5 dB at 32 MHz) against n = 6 for Butterworth.
Odd n also allows equal 50 ohm terminations.

**Ladder.** Matthaei closed-form g-values
(`g1 = 2a1/gamma`, `g_k = 4 a_{k-1} a_k / (b_{k-1} g_{k-1})`), denormalised as
`C = g/(R w_p)`, `L = g R / w_p`, shunt-C first: 679 pF, 1.22 uH, 1011 pF,
1.22 uH, 679 pF. Responses come from an ABCD cascade with
`S21 = 2 sqrt(Rs/RL) V_L/V_s`. For even n, `g_{n+1}` is a load conductance
when the last element is a series L; the n = 4 test catches that.

**Sallen-Key.** Unity-gain sections with equal R:
`w0 = 1/(R sqrt(C1 C2))`, `Q = sqrt(C1/C2)/2`, placed on the prototype poles,
plus a passive RC for the real pole. The Q = 4.55 section needs C1/C2 = 83.
Each stage is solved from its node equations with a single-pole op-amp
`A = wt/s`. To first order, finite GBW adds `2 Q w0 / wt` to the s^2
coefficient, which is why this section is so hungry.

**Dissipation.** With component Q, an element `x(s)` becomes `x(s + w/Q)`.
Perturbing ln H gives

```
delta ln|S21| = sum_k  (d phase / d ln x_k) / Q_k
```

and for uniform Q this collapses to `8.686 w tau_g / Q` dB, because scaling
every reactance is a frequency scaling. The rule predicts the inductor Q needed
for a dissipation budget at the band edge.

**Orchard.** A lossless ladder between matched terminations has |S21| <= 1,
and |S21| = 1 at each reflection zero `w = cos((2k-1) pi / 2n) w_p`. A maximum
cannot move to first order, so dS21/dx = 0 there for every element. The
active cascade has no such bound.

**Monte Carlo.** 20 000 trials per row, uniform tolerances, ripple = max - min
gain over 0-8 MHz, rejection = worst of 32, 36, 40, 60 MHz.

## Files

- `filters.py` - prototype, ladder and Sallen-Key models
- `run_filter.py` - design, dissipation, sensitivity, GBW and yield study
- `test_filters.py` - 21 checks

## Checks (`test_filters.py`)

- g-values against Matthaei Table 4.05-2 for 0.01/0.1/0.5 dB, n = 3, 4, 5, and
  Butterworth n = 5; order is minimal for both families; poles match `cheb1ap`.
- Ladder S21 equals `scipy.signal.cheby1` (complex, magnitude and phase) to
  1e-11 for n = 3, 5, 7 and for even n = 4 with its -0.5 dB DC level.
  Butterworth ladder is exactly -3 dB at the edge.
- |S11|^2 + |S21|^2 = 1 to 1e-12; S11 = 0 at the reflection zeros.
- Group delay from the pole formula matches the numerical derivative of the
  ladder phase to 1e-6.
- Uniform-Q loss matches `8.686 w tau / Q`, unequal Q matches the phase
  sensitivity sum, and the sensitivities add up to `-w tau`.
- Orchard: first-order sensitivity < 1e-8 at reflection zeros; 500 random
  +-10 % ladders never exceed |S21| = 1.
- Sallen-Key cascade equals the ladder to 1e-11; the finite-GBW stage matches
  an independent 3x3 nodal solve to 1e-10; the GBW search returns the 0.1 dB
  threshold.

## Run

```
pip install numpy scipy pytest
python3 run_filter.py
python3 -m pytest -q test_filters.py
```

## Sample output

```
== Ladder dissipation, QL = 60, QC = 1000 ==
  f MHz   delay ns  exact dB   predicted dB
   1.00     79.91    0.0266    0.0267
   4.00     89.71    0.1768    0.1778
   8.00    210.63    0.7804    0.8074
  QL for 0.25 dB extra at 8 MHz: 224.2 -> QL 225 gives 0.2467 dB exact

== First-order passband sensitivity, |dB| per 1 % component change ==
  7.608 MHz  ladder  4.82e-12   Sallen-Key  1.65e-01 (reflection zero)
  4.702 MHz  ladder  9.64e-12   Sallen-Key  6.29e-02 (reflection zero)
  3.200 MHz  ladder  8.59e-03   Sallen-Key  4.77e-02

== Op-amp for the cascade ==
  GBW 1.0 GHz: passband error 1.196 dB, 2 Q w0 / wt = 7.40 %
  GBW 4.0 GHz: passband error 0.306 dB, 2 Q w0 / wt = 1.85 %
  GBW for <= 0.1 dB passband error: 12306 MHz

== Monte Carlo, 20000 trials, uniform tolerances ==
  case                                   ripple p50    p99  rej min   yield
  LC ladder, L/C +-5 %                        0.655  1.679    72.45  84.78%
  LC ladder, L +-5 %, C +-2 %                 0.579  1.336    73.19  91.74%
  LC ladder, L +-2 %, C +-2 %                 0.557  0.887    73.65  99.98%
  LC ladder, L +-5 %, C +-2 %, QL 60          1.298  2.101    73.13  13.89%
  LC ladder, L/C +-2 %, QL 225                0.747  1.141    73.68  92.91%
  Sallen-Key, R +-1 %, C +-5 %                1.122  2.254    72.34  37.95%
  Sallen-Key, R +-1 %, C +-2 %                0.758  1.244    73.54  89.12%
  Sallen-Key, R +-1 %, C +-1 %                0.654  0.963    73.88  99.50%
  Sallen-Key, R/C +-1 %, 4 GHz GBW            0.798  1.215    74.06  86.25%
```

## Takeaways

- With the same 5 % capacitors the ladder yields 85 % against 38 % for the
  cascade. Its sensitivity vanishes at every ripple peak, so tolerance only
  pulls the passband down a little.
- The ladder's catch is dissipation. Ordinary chip inductors (Q around 60)
  add 0.8 dB at the band edge and drop yield to 14 %. Holding the edge loss to
  0.25 dB needs Q of about 225 (air-core or larger wirewound parts), or a
  predistorted design.
- The cascade needs tight capacitors and a very fast op-amp. A 0.1 dB budget
  needs about 12 GHz of GBW on the Q = 4.55 section, so in practice it would be
  predistorted for the chosen op-amp or built with a less GBW-sensitive
  topology.
- Rejection is not the issue for either: 72 dB or better in every trial,
  against a 60 dB requirement.
