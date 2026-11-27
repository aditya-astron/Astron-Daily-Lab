# 067 - Gravity-gradient boom sizing: stability regions and eccentric-orbit pitch libration

## Problem

A 50 kg microsat (bus inertia 2.2 / 2.5 / 1.6 kg m^2 roll/pitch/yaw) is held
nadir-pointing by gravity gradient alone, with a deployable boom and 1.5 kg tip
mass on the yaw axis. The orbit is 600 x 880 km, so e = 0.0197 and the local
vertical swings non-uniformly around the orbit, which drives pitch. How long
does the boom need to be to keep the pitch libration under 2 deg? Is the
result stable in all three axes, and does it stay clear of the parametric
resonances?

## Approach

**Reference model.** Full nonlinear rigid body in a Keplerian orbit:
quaternion attitude relative to LVLH, Euler's equations with the
gravity-gradient torque `3 mu/r^3 (c x I c)`, and true anomaly integrated
alongside. Every reduced model below is checked against it.

**Circular-orbit linear stability** (Hughes ratios `k1 = (I2-I3)/I1`,
`k2 = (I1-I3)/I2`, `k3 = (I2-I1)/I3`, time in 1/w0):

    pitch:     theta'' + 3 k2 theta = 0
    roll/yaw:  phi'' - (1-k1) psi' + 4 k1 phi = 0
               psi'' + (1-k3) phi' + k3 psi   = 0
               s^4 + (1 + 3k1 + k1 k3) s^2 + 4 k1 k3 = 0

This gives the Lagrange triangle `k1 > k3 > 0` (I2 > I1 > I3, a true minimum
of the Jacobi integral) and the thin DeBra-Delp lobe at negative k, which is
gyroscopically stable but sits on a saddle of the potential.

**Eccentric pitch** in true anomaly:

    (1 + e cos nu) theta'' - 2 e sin nu (theta' - 1) + 3 k2 sin theta cos theta = 0

- Forced response to first order: `theta = A sin nu`, `A = -2e/(3k2 - 1)`.
  It resonates at k2 = 1/3 and tends to e rad as k2 -> 1.
- Parametric stability comes from Floquet theory. The damping term has zero
  mean, so `det M = 1` and the motion is bounded iff `|tr M| < 2`. Removing the
  theta' term reduces the problem to a Mathieu equation with `delta = 3 k2`,
  `eps = e(1 - 3k2)`. The tongues sit at k2 = n^2/12, and the first one is
  `1/12 -+ e/8`.
- The exact 2pi-periodic libration is found by Newton shooting on the
  nonlinear equation and its variational equations. Its monodromy trace gives
  its stability.

**Sizing rule.** The boom must be long enough for 2 deg, plus a 15 % margin,
and the pitch frequency must stay at least 0.05 w0 away from 3/2 (the n = 3
tongue).

## Results

- **n = 1 tongue.** Edges match `1/12 -+ e/8` with an O(e^2) error.
- **n = 2 tongue.** It never opens. At k2 = 1/3 the monodromy is exactly the
  identity for every e (coexistence), so that point is a pure forced
  resonance.
- **n = 3 tongue.** The width grows as e^3: 1.1e-7 at e = 0.02 and 1.2e-4 at
  e = 0.2.
- **Amplitude limit alone.** Reaching 2 deg needs L = 1.34 m. That is
  w_p = 1.46 w0, and with a 15 % length margin L = 1.545 m, which lands almost
  exactly on the n = 3 tongue (L = 1.543 m, |tr M| = 1.997). The tongue is far
  too thin to cause growth by itself, but the design would have no detuning
  margin.
- **Design: L = 1.9 m.** k = (0.863, 0.806, 0.188), Lagrange region. Pitch is
  at 1.555 w0 and roll/yaw at 1.890 / 0.426 w0. The periodic libration is
  1.59 deg (linear 1.589) with tr M = -1.88. With +-10 % on the bus yaw
  inertia, w_p stays in 1.54-1.57.
- **Deployment.** The boom passes the k2 = 1/3 forced resonance at L = 0.19 m.
  The steady libration is already 30 deg at L = 0.25 m, so the deployment
  should not dwell near that length.
- **Full 3D model, 50 orbits.** Starting on the periodic libration, the pitch
  returns to 7e-14 rad. With 3 deg roll/yaw errors at the start, yaw peaks at
  9.2 deg because the yaw mode is soft (0.43 w0). Gravity gradient alone does
  not damp that, so a real design needs a damper.

## Verification (`cargo test --release`, 17 tests)

- **Kinematics.** The quaternion DCM matches the 3-2-1 rotation product, and
  quaternion rates reproduce `C' = -[w x] C` by finite difference.
- **Jacobi integral.** It is conserved to < 1e-9 over 10 orbits at 10 deg
  amplitude, with fourth-order convergence.
- **Linear roll/yaw and pitch against the nonlinear model.** Tested in the
  Lagrange and DeBra-Delp regions. The relative error falls as a^2, because
  the potential is even and has no quadratic term.
- **Analytic exponents.** The characteristic polynomial of the transition
  matrix `Phi(T)` (Faddeev-LeVerrier) equals `prod(lambda - exp(s_i T))` for
  stable, DeBra-Delp and unstable cases.
- **Reference frequencies.** A thin boom gives roll 2 w0 and yaw sqrt(k3) w0.
  Pitch is sqrt(3) w0 at k2 = 1.
- **Regions.** The Lagrange area is 1/8 of the k-square. DeBra-Delp stays
  bounded in the 3D model even though a pure roll offset lowers its
  potential. An unstable point diverges.
- **Planar eccentric equation.** It matches the 3D model to 1e-9 at e = 0.1
  and e = 0.3 (including at 0.4 rad pitch).
- **Floquet checks.** det M = 1, and at e = 0 tr M = 2 cos(2 pi sqrt(3 k2)).
  M(1/3, e) = I. The n = 1 edges have an O(e^2) error, and the n = 3 width
  ratio is 8 for 2x e. Growth inside the tongue equals ln(lambda) per orbit to
  1e-6.
- **Periodic libration.** Its peak equals |A| with an O(e^2) error, and it
  closes in the 3D model after 5 orbits.

## Run

    cargo run --release      # design study, see sample_output.txt
    cargo test --release

Pure Rust, no dependencies; the full study runs in about a second.
