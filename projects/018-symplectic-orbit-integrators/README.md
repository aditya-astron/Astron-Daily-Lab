# 018 - Symplectic vs. RK4 integration of long-arc Keplerian orbits

A Molniya orbit (a = 26562 km, e = 0.74, perigee 528 km altitude) is a hard
case for a fixed-step propagator: the satellite crosses perigee at 10 km/s and
spends most of its 12 h period crawling around apogee. Over a 50-day mission
arc, which integrator gives the best trajectory per force evaluation, and what
does each one get wrong?

Goal: implement RK4, leapfrog and Yoshida's 4th-order composition, compare them
against the exact two-body solution over 100 revolutions, and verify the
properties that distinguish them: order of accuracy, energy behaviour,
angular-momentum conservation and time reversibility.

## Approach

**Equations of motion.** `r'' = -mu r / |r|^3`, mu = 398600.4418 km^3/s^2.
Hamiltonian `H = |v|^2/2 - mu/|r|` splits into kinetic and potential parts,
each with an exact flow (drift and kick).

**Integrators.**

- RK4: classical 4-stage Runge-Kutta, 4 force evaluations per step. Accurate
  but not symplectic, so its energy error accumulates linearly with time.
- Leapfrog (kick-drift-kick): `v += h/2 a(r); r += h v; v += h/2 a(r)`.
  Composition of exact sub-flows, hence symplectic and time-reversible. It
  exactly conserves a "shadow" Hamiltonian `H + O(h^2)`, so the energy error is
  bounded forever. Kicks are radial and drifts are along v, so `r x v` is
  preserved to round-off. 1 force evaluation per step.
- Yoshida-4: leapfrog substeps of `w1 h, w0 h, w1 h` with
  `w1 = 1/(2 - 2^(1/3))`, `w0 = -2^(1/3)/(2 - 2^(1/3))`. The conditions
  `2 w1 + w0 = 1` and `2 w1^3 + w0^3 = 0` cancel the h^3 error term, giving a
  4th-order symplectic method at 3 force evaluations per step.

**Reference.** Exact propagation with Lagrange f and g coefficients in the
eccentric-anomaly difference dE, solving

    n dt = dE - (1 - r0/a) sin dE + (r0.v0 / sqrt(mu a)) (1 - cos dE)

by Newton iteration (whole revolutions removed first), then

    f = 1 - a/r0 (1 - cos dE)           g = dt - sqrt(a^3/mu) (dE - sin dE)
    fdot = -sqrt(mu a)/(r r0) sin dE     gdot = 1 - a/r (1 - cos dE)

## Verification (`test_integrators`)

- `from_elements` at perigee gives r = a(1-e), vis-viva speed, r.v = 0 and the
  correct inclination.
- Kepler propagator: returns to the start after one period (2e-11 km), runs
  37.3 revolutions forward and back (6e-8 km), matches an independent
  element-based position at M0 + n dt, conserves energy to 1e-13.
- Observed order from error halving on an e = 0.3 orbit: RK4 4.02,
  leapfrog 2.00, Yoshida-4 4.00 (tolerance 0.1). Yoshida weight conditions hold
  to machine precision.
- RK4's energy error sampled at the same orbital phase grows exactly linearly
  (rev 100 / rev 10 = 10.00).
- Leapfrog and Yoshida-4 energy envelopes are identical at rev 1 and rev 100.
- Angular momentum drift over 100 revolutions: under 1e-13 for both symplectic
  methods, 1.8e-8 for RK4.
- Leapfrog run forward 10 revolutions, velocity flipped, run forward again,
  lands back on the start within 2e-7 km; RK4 misses by 0.19 km.

## Run

    make
    ./test_integrators
    ./orbit_demo

## Sample output

```
method     N/rev   f-evals    max|dE/E|    max|dh/h|  final |dr| km
RK4         2000    800000    3.040e-07    1.817e-08     1.0394e+01
Leapfrog    2000    200000    4.700e-04    6.511e-14     2.0962e+04
Yoshida4    2000    600000    1.977e-07    9.472e-14     1.2795e+01

Energy error vs time, N = 2000 steps/rev (max |dE/E| within the revolution)
method           rev 1       rev 10      rev 100     100/10
RK4          1.706e-08    4.314e-08    3.040e-07       7.05
Leapfrog     4.700e-04    4.700e-04    4.700e-04       1.00
Yoshida4     1.977e-07    1.977e-07    1.977e-07       1.00

Equal cost, 12000 force evaluations per revolution, 100 revs
method     N/rev    max|dE/E|  final |dr| km
RK4         3000    4.106e-08     1.4353e+00
Leapfrog   12000    1.305e-05     8.0442e+02
Yoshida4    4000    1.237e-08     8.0041e-01
```

(The RK4 envelope ratio is 7 rather than 10 because the per-revolution maximum
includes a fixed periodic excursion at perigee on top of the linear drift.)

## Takeaways

- Bounded energy is not the same as an accurate trajectory. Leapfrog keeps
  energy and angular momentum perfectly bounded, yet after 100 revolutions it
  is 20 000 km off because the shadow orbit has a slightly wrong period, and
  that phase error grows linearly.
- At equal cost, Yoshida-4 beats RK4 slightly on position and by 3x on energy,
  and unlike RK4 it will not let the orbit decay or spiral out over very long
  arcs.
- For a highly eccentric orbit the real fix is not a better fixed-step method
  but variable steps or time regularisation (Sundman `dt = r ds`), which
  concentrates steps at perigee. Naive step-size adaptation breaks
  symplecticity; that trade-off is the next thing to look at.
