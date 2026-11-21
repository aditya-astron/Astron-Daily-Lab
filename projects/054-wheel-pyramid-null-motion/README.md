# 054 - Four-wheel pyramid: momentum envelope and zero-crossing avoidance

A bus with four reaction wheels in a pyramid has one more wheel than it needs.
Two questions follow. How much momentum can the cluster store in each
direction, and at what pyramid elevation is the worst direction best? And how
should the spare degree of freedom be used? The standard minimum-norm
distribution parks the wheels near zero speed, where bearing stiction,
lubricant starvation and tachometer sign ambiguity live. Spending the spare
freedom on a speed bias avoids that, and the question is what it costs.

## Approach

**Geometry.** Wheel axis `a_i = [cos b cos phi_i, cos b sin phi_i, sin b]`,
`phi_i = 45, 135, 225, 315 deg`, body momentum `h = A x`, `|x_i| <= hmax`.
The reachable set is the zonotope `Z = A [-hmax, hmax]^4`.

**Envelope radius, two independent ways.**

- LP: `max t  s.t.  A x = t d,  |x_i| <= hmax` (HiGHS).
- Facets: a 3D zonotope's facet normals are `+-(a_i x a_j)`, its support
  function is `hZ(n) = hmax sum_k |a_k . n|`, so
  `r(d) = min_{n.d > 0} hZ(n) / (n.d)`.

The inscribed sphere, the capacity guaranteed in every direction, is the
distance to the nearest facet plane: `min_n hZ(n) / |n|`.

**Closed forms used as checks.**

- Along z all wheels contribute `sin b`: `r_z = 4 sin b`.
- Along x, `x = [1, -1, -1, 1] hmax` gives `r_x = 2 sqrt2 cos b`; the null
  vector can only raise its largest entry.
- At `tan b = 1/sqrt2` (35.26 deg) the axes are cube diagonals
  `(+-1, +-1, 1)/sqrt3`, `A A^T = 4/3 I`. The facet spanned by `(1,1,1)` and
  `(-1,1,1)` has normal `(0,-1,1)/sqrt2`; the other two wheels project
  `2/sqrt6` on it, so the inscribed sphere is `4/sqrt6 hmax = 1.633 hmax`.
  A bounded search over b lands on this angle.
- Minimum-norm law there: `A^+ = 3/4 A^T`, so along wheel 1
  `A^+ a1 = [3/4, 1/4, -1/4, 1/4]` and it reaches `4/3 hmax`, while
  `x = [1, 1, -1, 1]` reaches `2 hmax` (`a2 - a3 + a4 = a1`). The pseudoinverse
  alone gives away a third of the envelope there, and nowhere more.

**Null-motion bias.** `A n = 0` with `n = [1, -1, 1, -1]` for every b, so
`w = A^+ h + lam n` leaves the body untouched for any `lam`. Since
`n_i w_i = n_i p_i + lam`, all wheels sit in `[hmin, hmax]` (with signs `n`)
when `lam` is in `[hmin - min(n p), hmax - max(n p)]`; the midpoint centres
the margins. The set of `h` this can hold is `A diag(n) [hmin, hmax]^4`, whose
centre `A n (hmin+hmax)/2` is zero. It is therefore the full envelope scaled
by `(hmax - hmin) / (2 hmax)`: never crossing zero costs a bit more than half
the capacity.

**Run.** 4 N m s / 6000 rpm wheels, 600 rpm floor, three LEO orbits of an
inertially held bus with orbit-rate and twice-orbit-rate gravity-gradient
momentum. `lam` tracks the target through a 1 Hz saturated first-order law
with 0.02 N m per wheel reserved for redistribution, starting from all
wheels at rest.

## Verification

`test_wheel_pyramid.py` (18 tests):

- `A n = 0` at every elevation.
- LP radius equals `4 sin b` and `2 sqrt2 cos b` on the axes to 1e-9.
- LP equals the facet formula in 200 random directions per elevation.
- Facet-distance sphere equals the minimum over 6000 sampled directions.
- Optimum elevation equals `atan(1/sqrt2)` and the sphere `4/sqrt6`.
- Minimum-norm reach `4/3` versus `2` along a wheel axis, and the 2/3
  ratio is the worst over 5000 directions.
- Bias-constrained LP equals the scaled envelope in 100 directions.
- Bias target puts every wheel inside the band with equal margins.
- Orbit run: zero crossings without the bias, none with it after spin-up.
  `|A w - h| < 1e-12` throughout, and the null torque stays within its
  allowance.

## Run

```
pip install numpy scipy pytest
python3 run_study.py
python3 -m pytest -q
```

## Sample output

```
Envelope vs pyramid elevation (units of hmax)
 beta[deg]  r_x(LP)  2r2*cos  r_z(LP)  4*sin   sphere
   20.000   2.6579   2.6579   1.3681  1.3681  1.2164
   30.000   2.4495   2.4495   2.0000  2.0000  1.5492
   35.264   2.3094   2.3094   2.3094  2.3094  1.6330
   45.000   2.0000   2.0000   2.8284  2.8284  1.4142
   60.000   1.4142   1.4142   3.4641  3.4641  1.0000

best elevation 35.2644 deg (atan(1/sqrt2) = 35.2644), sphere 1.632993 hmax (4/sqrt6 = 1.632993)
along wheel 1 axis: envelope 2.0000, min-norm 1.3333 hmax (closed forms 2 and 4/3)
envelope over 5000 directions: 1.6330 .. 2.2822 hmax
min-norm law alone reaches 0.6700 of the envelope at d = [-0.572 -0.585 -0.575]
zero-crossing-free envelope = 0.450 x full envelope; guaranteed 2.939 N m s in any direction

Three-orbit run, wheels 4.0 N m s at 6000 rpm, floor 600 rpm, null torque 0.02 N m
  min-norm only     zero crossings  29  min |w|    0.0 rpm  max |w| 2329.4 rpm  worst wheel below floor  49.8 %  |Aw-h| 6.7e-16
  null-motion bias  zero crossings   0  min |w| 1630.9 rpm  max |w| 5068.0 rpm  worst wheel below floor   0.0 %  |Aw-h| 1.1e-15
```

The minimum-norm law spends half the time with at least one wheel below
600 rpm and crosses zero 29 times in three orbits. The bias keeps every wheel
above 1600 rpm at the cost of running them faster, which is the usual trade
of bearing life and power against stiction and control glitches.

The minimum-norm worst direction is sampled close to `-a1`. The sampled
envelope maximum is 2.28 rather than the axis value 2.31 because no sample
falls exactly on an axis.
