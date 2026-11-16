# 040 - Earth-Moon L2 Halo Family to the 9:2 NRHO

**Domain:** Orbital mechanics

## Problem

Gateway-class missions fly a near-rectilinear halo orbit (NRHO) about
Earth-Moon L2 whose period is in 9:2 resonance with the synodic month, so
the orbit never falls into Earth eclipse for long. That orbit is a member
of the L2 southern halo family, which itself branches off the planar
Lyapunov family. The goal is to build that whole chain in the circular
restricted three-body problem (CR3BP): find the branch point, seed a halo
from analytic theory, continue the family through its fold, and land on the
orbit with T = 29.530589 / 4.5 = 6.5624 days, then report its perilune,
apolune and stability.

## Approach

**Model.** Barycentric rotating frame, mu = GM_moon / (GM_earth + GM_moon)
= 0.0121506, unit length 384,400 km, unit time 4.3425 d.

    x'' - 2y' = U_x,  y'' + 2x' = U_y,  z'' = U_z
    U = (x^2 + y^2)/2 + (1-mu)/r1 + mu/r2,     C = 2U - v^2

L1 and L2 come from the collinear quintics in gamma. The 6x6 state
transition matrix is integrated alongside the state (DOP853, tol 1e-12).

**Richardson third-order seed** (`richardson.py`). Lindstedt-Poincare
expansion about L2 in gamma-scaled coordinates. Removing secular terms at
third order gives the amplitude constraint

    l1 Ax^2 + l2 Az^2 + Delta = 0,   Delta = lambda^2 - c2

so halos only exist above Ax_min = sqrt(-Delta/l1).

**Differential correction** (`families.py`). Every orbit is symmetric about
the x-z plane, so it starts at a perpendicular crossing
(x0, 0, z0, 0, vy0, 0) and is periodic when the next crossing has
vx = vz = 0. The crossing time floats, so the Jacobian is

    d s(T/2) / dX = Phi[:, X] - f(s) Phi[1, X] / vy

**Pseudo-arclength continuation.** The family tangent n is the null vector
of the 2x3 constraint Jacobian in (x0, z0, vy0). Each step adds
(X - X*) . n = ds as a third equation, which carries the solver round the
fold where Az peaks (~77,600 km) and starts falling again towards the NRHO.
Once the period drops below 6.5624 d the last orbit is re-solved with
T/2 = T_target as the third equation.

**Monodromy and stability.** For symmetric orbits
M = G Phi(T/2)^-1 G Phi(T/2) with G = diag(1,-1,1,-1,1,-1). The inverse
uses the symplectic identity Phi^-1 = Omega^-1 Phi^T Omega instead of a
numerical inverse. Eigenvalues come in reciprocal pairs plus a trivial pair
at +1; each nontrivial pair gives a stability index nu = (lam + 1/lam)/2.

**Bifurcation.** On a planar Lyapunov orbit the (z, vz) block of M is
decoupled with det 1, so the halo branch appears where its trace equals 2.
That point is found with Brent's method along the Lyapunov family.

## Verification (`python3 test_halo_nrho.py`, 12 checks)

| check | result |
|---|---|
| L1, L2, C(L1), C(L2) vs. reference 0.836915, 1.155682, 3.18834, 3.17216 | < 2e-6 |
| Richardson Sun-Earth L1 coefficients (c2..c4, lambda, k, Delta, s1, s2, l1, l2) vs. his 1980 table | < 2e-6 rel |
| ISEE-3 class halo: Az 120,000 km -> Ax ~ 206,000 km, T ~ 177.7 d | yes |
| Richardson velocity = d/dt of Richardson position | 1e-9 |
| STM vs. central finite differences | < 1e-6 rel |
| STM symplectic, det = 1 | 1e-9 |
| half-period monodromy vs. full-period propagated STM | < 1e-6 rel |
| bifurcation: C = 3.1521, T = 14.83 d, Ax within 2 % of Richardson Ax_min | yes |
| halo C and T approach the bifurcation quadratically in Az (ratio 4 for 2x Az) | yes |
| 9:2 NRHO: rp 3,249 km, ra 71,222 km, nu = -1.32 (published ~3,200-3,500, ~71,000, ~1.3) | yes |
| reciprocal eigenvalue pairs, trivial pair at +1 | yes |
| Jacobi constant over 3 NRHO periods | < 1e-10 |
| every continuation step lands on vx = vz = 0 with step length ds | yes |

## Run

    pip install numpy scipy
    python3 halo_nrho.py          # ~11 s
    python3 test_halo_nrho.py     # or: python3 -m pytest test_halo_nrho.py

## Sample output

    Lyapunov -> halo bifurcation (out-of-plane pair at +1)
      C = 3.152119, T = 14.8319 d, Ax = 11630 km (Richardson Ax_min 11430 km)

    Richardson guess at Az = 10000 km, corrected at fixed z0
      dx0 = 493.3 km, dvy0 = 0.59 m/s, T 14.7753 -> 14.7997 d

    Pseudo-arclength continuation (37 orbits)
       z0 [km]   rp [km]   ra [km]   T [d]      C        nu1      nu2
        -11557     50330     74757  14.800   3.14819  566.055   0.985
        -45220     42802     81011  14.278   3.09840  207.904   0.506
        -74976     25303     88182  11.958   3.02292   10.841  -1.000
        -77583     14317     84066   9.635   3.01635   -0.802  -0.079
        -74676      7500     77802   7.889   3.02829   -1.692   0.364
        -70000      3249     71222   6.562   3.04649   -1.323   0.683

    9:2 NRHO (southern, perilune over the north pole)
      s0 = [1.0220282130, 0.0000000000, -0.1821013936, 0.0000000000, -0.1032709456, 0.0000000000]
      T = 1.51119943 = 6.56235 d,  C = 3.04649375
      perilune radius 3249 km (altitude 1512 km), apolune radius 71222 km
      stability indices -1.3230, 0.6829;  max |lambda| 2.1892
      one-period closure error 1.1 mm

The family is wildly unstable near the bifurcation (nu ~ 566), passes
through a linearly stable window around rp ~ 14,000 km, and the 9:2 NRHO
sits just outside it with a growth factor of 2.2 per 6.6-day revolution,
which is why it needs only small, infrequent station-keeping burns.

## Files

- `cr3bp.py` - dynamics, STM, libration points, crossing events, monodromy
- `richardson.py` - third-order halo approximation
- `families.py` - single shooting, pseudo-arclength continuation, stability
- `halo_nrho.py` - the full pipeline and report
- `test_halo_nrho.py` - verification
