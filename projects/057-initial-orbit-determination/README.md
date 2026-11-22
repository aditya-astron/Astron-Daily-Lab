# 057 - Initial orbit determination: Gibbs, Herrick-Gibbs and Gauss angles-only

A new object shows up on a radar or in a telescope frame and there is no
orbit to start a filter or batch fit from (project 026 assumed one). Initial
orbit determination gets a full state from three sightings and nothing else.
This project implements the three classical methods, checks them against
textbook examples and exact synthetic data, then measures the things the
textbooks only hint at: where Gibbs and Herrick-Gibbs swap places under
noise, how fast angles-only accuracy collapses as the sightings get closer,
and where the standard Gauss refinement stops converging.

## Approach

**Propagation.** Exact two-body motion by universal variables: Newton on
the universal Kepler equation (its derivative is the radius), Stumpff
functions with series near z = 0, Lagrange f, g, fdot, gdot. Checked
against DOP853 at rtol 1e-13 for near-circular, eccentric, hyperbolic and
backwards cases, plus the identity f gdot - fdot g = 1.

**Gibbs** (three positions, no times):

```
N = r1 (r2 x r3) + r2 (r3 x r1) + r3 (r1 x r2)
D = r2 x r3 + r3 x r1 + r1 x r2
S = r1 (r2 - r3) + r2 (r3 - r1) + r3 (r1 - r2)
v2 = sqrt(mu / (|N| |D|)) (D x r2 / r2 + S)
```

**Herrick-Gibbs** (three positions with times): Taylor expansion of r(t)
about t2 with the two-body acceleration substituted,

```
v2 = -dt32 (1/(dt21 dt31) + mu/(12 r1^3)) r1
     + (dt32 - dt21)(1/(dt21 dt32) + mu/(12 r2^3)) r2
     + dt21 (1/(dt32 dt31) + mu/(12 r3^3)) r3
```

**Gauss angles-only.** Sightlines L1..L3 from site vectors R1..R3 (oblate
Earth). With f and g truncated after the mu/r^3 term, the middle radius
solves

```
x^8 + a x^6 + b x^3 + c = 0,   a = -(A^2 + 2AE + R2^2),  b = -2 mu B (A + E),  c = -(mu B)^2
```

where A, B come from the D_ij = R_i . (L_j x L_k) determinants and
E = R2 . L2. Roots below the surface or with negative range are rejected.
The ranges then follow linearly, and v2 from r1, r3 and the series f, g.
Two refinements:

- *Fixed point* (Curtis algorithm 5.6): recompute exact f, g from the
  current (r2, v2), average with the previous ones, re-solve the ranges.
- *Newton correction*: three angle pairs fix the six-element state
  exactly, so solve for (r2, v2) that reproduces all three sightlines,
  with a central-difference Jacobian and step halving. Same answer where
  the fixed point converges, but it keeps converging on long arcs.

## Results

**Reference examples.** Gibbs reproduces Curtis 5.1 to all printed digits;
Herrick-Gibbs matches Vallado 7-4 to 0.1 m/s. Gauss on Curtis 5.11/5.12
lands 1.2 km and 2 m/s from the book. Rounding the example's inputs by
half a printed digit moves r2 by 1.0 km median (2.8 km max), so the gap is
the data, not the method; on exact synthetic data Gauss recovers the orbit
to under a millimetre.

**Gibbs vs Herrick-Gibbs** (700 km SSO, 10 m per-axis position noise):

| separation | HG truncation | Gibbs RMS | HG RMS |
|---|---|---|---|
| 1 deg  | 1.6e-5 m/s | 41.7 m/s | 0.75 m/s |
| 5 deg  | 9.9e-3 m/s | 1.70 m/s | 0.15 m/s |
| 12 deg | 0.33 m/s   | 0.31 m/s | 0.33 m/s |
| 30 deg | 12.3 m/s   | 0.06 m/s | 12.3 m/s |

- Herrick-Gibbs truncation is fourth order in the spacing (measured order
  4.00 +/- 0.15, tested). Its noise response sits on the central-difference
  floor sigma sqrt(3/2)/dt.
- Gibbs is exact on noise-free data at any separation, but its cross
  products of nearly parallel vectors amplify noise roughly as 1/sep^2.
- They swap at about 12 deg here, further out than the 1-5 deg rule of
  thumb, because with 10 m radar noise the HG truncation stays below the
  noise for longer.

**Gauss on a LEO pass** (800 km, 65 deg, 520 s pass from a 24.65 N,
1.68 km site, 1 arcsec noise, Newton refinement, medians of 500 trials):

| spacing | series start error | |dr2| | |dv2| | |da| |
|---|---|---|---|---|
| 10 s  | 0.03 km | 13.9 km | 109 m/s | 233 km |
| 40 s  | 0.50 km | 0.92 km | 7.3 m/s | 15.7 km |
| 80 s  | 2.0 km  | 0.20 km | 1.6 m/s | 3.3 km |
| 160 s | 8.0 km  | 0.06 km | 0.55 m/s | 1.2 km |

The range comes from the curvature of the track, so position error falls
as 1/spacing^2; spread the sightings across the pass. The series start
error grows as spacing^2 the other way, which is what refinement removes.

**Gauss on a GEO object** (i = 3 deg, 1 arcsec): the fixed-point iteration
needs 76-158 iterations up to 1 h spacing and diverges from 90 min on;
the Newton correction converges in 3-6 iterations at every spacing out to
4 h (120 deg of arc), with a 0.4 km median position error there. At 6 h
spacing (180 deg of arc) the series start has no physical root at all.

## Run

```
pip install numpy scipy pytest
python3 run_iod.py      # ~1 min, prints sample_output.txt
pytest -q               # 15 tests, ~5 s
```

## Files

- `iod.py`: propagator, elements, Gibbs, Herrick-Gibbs, site geometry, Gauss, Newton angles-only correction
- `run_iod.py`: the study
- `test_iod.py`: textbook, ODE, convergence-order and recovery checks
- `sample_output.txt`: output of `run_iod.py`

Two-body only: no J2, light time, aberration or refraction, which matter
for a real catalogue but not for comparing the methods.
