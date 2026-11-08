# 022 - Minimum-Length Planar Nozzle by the Method of Characteristics

**Domain:** Rocket propulsion

## Problem

Design the shortest supersonic nozzle that turns a sonic throat flow into
uniform, perfectly axial flow at a chosen exit Mach number, with no shocks.
Then check the contour three independent ways and compare it with the simple
straight-walled alternative a designer would otherwise pick.

## Approach

Planar, steady, inviscid, homentropic, irrotational flow. The throat is a
sharp corner (the limiting case that gives minimum length). A centred
Prandtl-Meyer fan at the corner turns the wall to

```
theta_max = nu(Me) / 2
nu(M) = sqrt((g+1)/(g-1)) atan( sqrt((g-1)/(g+1) (M^2-1)) ) - atan( sqrt(M^2-1) )
```

The fan waves reflect off the centreline and the wall downstream is shaped to
absorb each one, so nothing reflects back into the core. Along the
characteristics the Riemann invariants are constant:

```
C- (slope theta - mu):  K- = theta + nu
C+ (slope theta + mu):  K+ = theta - nu        mu = asin(1/M)
```

With `n` fan waves at `theta_j = nu_j`, the state at the crossing of reflected
C+ line `i` and fan line `j` is algebraic: `theta = theta_j - theta_i`,
`nu = theta_j + theta_i`. Only the geometry needs the unit processes: each node
is the intersection of the C+ and C- lines from its two upstream neighbours,
using the average of the characteristic angles at both ends. Wall points sit
where each C+ line reaches the previous wall segment, extended at the mean of
the two wall-point flow angles.

**Verification (`nozzle_test.go`):**

1. `nu(M)`, `A/A*`, `p/p0` against NACA 1135 tables; inverse `nu -> M` to 1e-9.
2. Mass: exit half-height converges to the isentropic `A/A*` as `n` grows
   (error more than halves per doubling; 0.017% at n = 240).
3. Exit state: wall and centreline reach exactly `Me` at `theta = 0`; the wall
   is monotone with no inflection; every net point lies inside the nozzle.
4. Axial momentum: throat momentum `p*(1+g)h*` plus the wall pressure integral
   `int p dy` must equal the exit flux `pe(1+g Me^2)he`. This closes only if
   both the contour and its wall pressures are right; residual falls with `n`.
5. Thrust coefficient: the momentum form equals the closed-form ideal Cf
   (Sutton eq. 3-30) to 1e-12; the MOC wall integral gives Cf_vac = 1.708 for
   gamma = 1.22, Me = 3.2, matching a hand calculation.
6. Planar wedge divergence factor `sin(a)/a` = 0.98862 at 15 deg.

## How to run

Requires Go 1.24+, no external modules.

```
go test -v ./...
go run ./cmd/nozzle
go run ./cmd/nozzle -csv wall.csv     # gamma = 1.22 contour: x, y, theta, M
```

## Sample output

```
Case 1: Me = 2.4, gamma = 1.40
  nu(Me) = 36.7465 deg, theta_max = 18.3733 deg, isentropic A/A* = 2.40310

     n    h_e/h*   area err    L/h*   momentum resid
     7   2.36348    -1.649%   7.9904   +8.71e-03
    15   2.38760    -0.645%   8.0503   +3.33e-03
    30   2.39671    -0.266%   8.0722   +1.35e-03
    60   2.40051    -0.108%   8.0812   +5.45e-04
   120   2.40206    -0.043%   8.0849   +2.18e-04
   240   2.40268    -0.017%   8.0863   +8.68e-05
   480   2.40293    -0.007%   8.0869   +3.45e-05

  observed order 1.33, extrapolated L/h* = 8.0873

Case 2: Me = 3.2, gamma = 1.22, n = 200
  theta_max      = 33.575 deg
  eps (A/A*)     = 8.3041 isentropic, 8.3002 MOC
  pe/p0          = 0.01524
  length L/h*    = 33.794
  Cf vac         = 1.7074 (throat + wall pressure), 1.7079 ideal 1-D

  15 deg planar wedge, same eps:
    length L/h*  = 27.259 (81% of MOC)
    lambda       = 0.9886
    Cf vac       = 1.6899 (1.05% below MOC)
```

## Notes

- Case 1 is Anderson's Example 11.1. The area error converges at roughly order
  1.3 because the first fan wave starts a full step off sonic; the contour
  length settles at about 8.09 throat half-heights.
- The 15 deg wedge is shorter here, not longer. At Me = 3.2 and gamma = 1.22
  the corner turns 33.6 deg and the wall then needs a long run to straighten
  the flow back to axial. "Minimum length" is minimum among shock-free nozzles
  with uniform exit flow; the wedge buys its length with a 1% Cf loss from
  divergence. Real engines truncate or Rao-optimise the MOC contour to sit
  between the two.
- Planar only. Axisymmetric MOC adds a source term so the invariants are no
  longer constant along characteristics, and needs a smooth throat; the
  planar case is the clean one to validate against.
