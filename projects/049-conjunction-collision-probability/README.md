# 049 - Conjunction collision probability and avoidance burn sizing

A 700 km sun-synchronous satellite gets a conjunction warning against a
tracked debris fragment: 62 m predicted miss, 72 deg crossing angle,
8.8 km/s relative speed, along-track uncertainties of hundreds of metres.
The questions are what the probability of collision really is, whether that
number can be trusted, and how small an along-track burn brings it below the
1e-4 manoeuvre threshold for a given lead time.

## Approach

**Encounter plane.** At 8.8 km/s the error ellipsoid passes in about 0.3 s,
so relative motion is a straight line and the position error is frozen
(short-term encounter). The frame has z along the relative velocity, x along
the miss vector; projecting the combined ECI covariance `C1 + C2` onto x-y
leaves a 2D Gaussian `N(m, C)` and

```
Pc = integral over |x| <= R of N(x; m, C) dA,    R = 15 m hard-body radius
```

**Four evaluations of the same integral.**

- Foster: polar quadrature, Gauss-Legendre in r and trapezoid in theta. The
  theta integrand is periodic and smooth, so convergence is exponential.
- Alfano 1D: in the principal frame,
  `Pc = int phi((x-mx)/sx)/sx [Phi((h-my)/sy) - Phi((-h-my)/sy)] dx`,
  `h = sqrt(R^2 - x^2)`, with `x = R sin t` to remove the chord-end root. The
  CDF difference is taken in whichever tail it lives in; done naively it
  cancels to 1 - 1 for far-off chords and loses every digit at Pc ~ 1e-15.
- Chan: `Pc = e^{-v/2} sum (v/2)^m/m! [1 - e^{-u/2} sum_{k<=m} (u/2)^k/k!]`,
  `u = R^2/(sx sy)`, `v = (mx/sx)^2 + (my/sy)^2`. This is the noncentral
  chi-square CDF; exact for isotropic C, an equal-area approximation otherwise.
- Small object: `pi R^2 N(0; m, C)`.

**Monte Carlo.** Samples relative position *and* velocity in ECI (velocity
error 0.15 m/s per axis per object), finds each sample's own straight-line
closest approach and counts misses under R. It never uses the projection, so
it checks the frame construction as well as the integral.

**Dilution.** Scaling `C -> k^2 C`, the small-object form
`Pc(k) = R^2/(2 k^2 sqrt|C|) exp(-d^2/(2k^2))` peaks at `k^2 = d^2/2` with
`Pc_max = R^2/(e d^2 sqrt|C|)`, d the Mahalanobis miss. If the real covariance
is larger than the peak scale, a low Pc can mean "we don't know where it is"
rather than "it will miss".

**Avoidance.** An impulse `dv` in RTN at lead time t moves the primary by the
Clohessy-Wiltshire solution

```
x =  sin(nt)/n vR + 2(1-cos nt)/n vT
y = -2(1-cos nt)/n vR + (4 sin(nt)/n - 3t) vT
```

and the new miss is the encounter-plane projection of `dr - d`. The smallest
along-track `dv` (either sign) reaching Pc = 1e-4 is bracketed on a log grid
and refined with Brent. Each answer is re-checked by integrating the burn in
two-body dynamics (DOP853, both arcs from the same back-propagated state so
integrator error cancels) and recomputing Pc from that displacement.

## Results

- All four methods agree with the Monte Carlo (z = -0.31 on 4e6 samples).
  Chan is 0.2 % high at aspect ratio 5.7; small-object is 0.4 % high since
  R/sigma_min = 0.23.
- Curvature over the encounter is 1e-6 m against a 15 m body: the
  rectilinear assumption is not the limiting error.
- Pc = 3.7e-3 sits on the diluted side (peak 8.2e-3 at k = 0.44). Shrinking
  the covariance with better tracking would roughly double Pc before it fell.
- The narrow axis of the projected ellipse is exactly the radial direction
  (sigma = sqrt(25^2 + 60^2) = 65 m). A burn half an orbit ahead converts
  along-track dv into radial offset `4 dv / n`, measured against 65 m instead
  of 371 m, so 35 mm/s at half an orbit beats 66 mm/s at a full orbit. At whole
  orbits the shift is purely along-track (`-3 dv t`) and dv falls as 1/lead.
- CW and two-body displacements differ by at most 0.11 m out of 1.2 km, and
  the two-body Pc after the burn stays below 1e-4.

## Run

```
pip install numpy scipy pytest
python3 run_conjunction.py
python3 -m pytest -q          # 32 checks, about 1 s
```

## Files

- `conjunction.py` - frames, four Pc methods, Monte Carlo, scaling, CW and two-body burn models
- `run_conjunction.py` - the scenario study
- `test_conjunction.py` - checks

## Checks (`test_conjunction.py`)

- Hand value: isotropic centred sigma = 100 m, R = 20 m gives
  `1 - e^{-0.02} = 0.0198013266932` from all three integrators.
- Isotropic offset cases against `scipy.stats.ncx2.cdf` to 1e-9, including
  tails down to 2e-11 and R larger than sigma.
- Foster vs. Alfano on 12 random ellipses up to aspect 30 and Pc down to
  3e-15 (the case that exposed the CDF cancellation), to 1e-9.
- Foster converges exponentially (1e-14 at 16 x 32 nodes for R = 2.5 sigma_min).
- Rotation invariance; Chan exact at aspect 1 and worsening with aspect;
  small-object ratio approaches 1 monotonically as R -> 0; Pc -> 1 for large R.
- Numerical maximum over k matches the closed-form small-object peak to 1e-4.
- Encounter frame orthonormal and right-handed, projected miss equals the
  rectilinear closest approach, RTN covariance rotation preserves the trace.
- Scenario Monte Carlo within 4 sigma of the analytic Pc.
- CW at whole orbits is pure `-3 dv t` drift; CW vs two-body for radial,
  along-track and cross-track burns within the `|d|^2/a` linearisation error.

## Sample output

```
Conjunction at TCA
  altitude                 700 km, crossing angle 72 deg
  relative speed           8821.8 m/s
  miss distance            62.18 m (rectilinear check 62.18 m, t = -1.2e-14 s)
  projected sigmas         370.7 m x 65.0 m (aspect 5.70)
  Mahalanobis miss         0.658
  hard-body radius         15.0 m

Pc by method
  Foster polar quadrature  3.744810e-03
  Alfano 1D integral       3.744810e-03   rel diff 4.6e-16
  Chan equivalent-area     3.753239e-03   rel err +2.25e-03
  small-object             3.760108e-03   rel err +4.09e-03
  Monte Carlo 3D (n=4,000,000)  3.735250e-03   14941 hits, z = -0.31

Short-term encounter assumption
  +/-5 sigma transit time  307.5 ms
  relative accel           1.08e-04 m/s^2
  curvature over transit   1.3e-06 m  (vs HBR 15 m)

Covariance scaling C -> k^2 C
  max Pc numeric           8.1760e-03 at k = 0.4435
  small-object closed form 7.9363e-03 at k = 0.4652
  k      Pc
  0.250  3.804e-03
  0.444  8.176e-03
  0.500  7.984e-03
  1.000  3.745e-03
  2.000  1.104e-03
  4.000  2.877e-04
  8.000  7.269e-05
  actual covariance (k=1) is above the peak scale: Pc is diluted; a lower Pc here is not reassurance

Avoidance: smallest along-track burn for Pc < 1e-04
  lead [orbits]  dv [mm/s]   along-track shift [m]  CW vs two-body [m]  Pc after
            0.5     -34.86                 +309.9               0.005  9.998e-05
            1.0     -66.28                +1178.4               0.107  9.993e-05
            1.5     -26.93                 +718.1               0.036  9.988e-05
            2.0     -33.14                +1178.4               0.100  9.992e-05
            5.0     -13.26                +1178.4               0.098  9.991e-05
```
