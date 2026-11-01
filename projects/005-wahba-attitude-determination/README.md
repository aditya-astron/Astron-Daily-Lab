# 005 - Wahba's Problem: TRIAD, Davenport q-method and QUEST

**Domain:** Attitude determination

## Problem

A spacecraft measures a few unit vectors in its own body frame (star
directions, the Sun, the geomagnetic field) and knows the same vectors in an
inertial reference frame from a catalogue or model. Find the rotation between
the two frames, and know how good that answer is.

Three classical solvers are implemented from scratch in dependency-free Rust
and compared on two sensor suites:

1. A star tracker: six stars in an 8 deg half-angle field of view, 3 arcsec
   noise per star.
2. A coarse suite: fine sun sensor (0.05 deg) and magnetometer (1 deg),
   60 deg apart.

## Approach

Conventions: `b = A r`, scalar-last quaternion,
`A(q) = (q4^2 - |qv|^2) I + 2 qv qv^T - 2 q4 [qv x]`, attitude errors in the
body frame, `A_est A_true^T = I - [dtheta x]`.

Wahba loss, weights `a_i ~ 1/sigma_i^2` normalised to sum to 1:

```
L(A) = 1/2 sum a_i |b_i - A r_i|^2 = 1 - tr(A B^T),   B = sum a_i b_i r_i^T
```

- **TRIAD.** Build an orthonormal triad from (v1, v1 x v2) in each frame,
  `A = [tb][tr]^T`. Trusts the first vector fully, so order matters.
- **Davenport q-method.** With `S = B + B^T`, `sigma = tr B`,
  `z = [B23-B32, B31-B13, B12-B21]`,
  `K = [[S - sigma I, z], [z^T, sigma]]`; the optimal quaternion is the
  eigenvector of K for its largest eigenvalue. Solved with cyclic Jacobi.
- **QUEST.** Newton on the characteristic polynomial of K starting from
  `lambda0 = sum a_i = 1` (converges in 1-2 steps), then
  `X = (alpha I + beta S + S^2) z`, `gamma = (lambda + sigma) alpha - det S`,
  `alpha = lambda^2 - sigma^2 + tr adj S`, `beta = lambda - sigma`,
  `q ~ [X, gamma]`.
  This is singular at 180 deg. Rotating the references 180 deg about axis i
  turns `tr B` into `K_ii`, so the frame to solve in is picked up front as
  the largest diagonal entry of K, with a normalised-gamma check as backup.

Covariances checked against the Monte Carlo:

```
QUEST / q-method:  P = [ sum sigma_i^-2 (I - b_i b_i^T) ]^-1
TRIAD (anchor 1):  P = s1^2 I + |b1 x b2|^-2 [ s1^2 (b1.b2)(b1 b2^T + b2 b1^T)
                                               + (s2^2 - s1^2) b1 b1^T ]
```

## Verification (`cargo test`)

- Quaternion -> DCM matches a hand-written 90 deg rotation; round trip over
  2000 random attitudes plus all Shepperd branches.
- Jacobi recovers the spectrum of `Q diag(4, -1, 2.5, 0.3) Q^T`.
- Noise-free data: all three solvers recover the attitude to < 1e-12 rad.
- QUEST and q-method agree to < 1e-9 rad and on `lambda_max` to 1e-12 over
  2000 noisy cases; minimum loss equals `1 - lambda_max`; perturbing the
  q-method answer always raises the loss.
- QUEST correct at exactly 180 deg and just short of it about five axes.
- Covariance hand case: two orthogonal vectors, equal sigma, gives
  `diag(s^2, s^2, s^2/2)`.
- Monte Carlo (20 000 trials) 1-sigma per axis within 3 % of both analytic
  covariances.

## Run

```
cargo run --release
cargo test --release
```

## Sample output

```
Case 1: star tracker, 6 stars, 8 deg half-FOV, 3 arcsec/star, 20000 trials
  1-sigma error, arcsec                 x         y   z(roll)
  predicted (QUEST covariance)     1.267     1.265    13.741
  QUEST                            1.267     1.265    13.840
  q-method                         1.267     1.265    13.840
  QUEST vs q-method max difference: 1.57e-6 arcsec, Newton iterations <= 2

Case 2: sun sensor 0.05 deg + magnetometer 1 deg, 60 deg apart, 20000 trials
  1-sigma error, deg                    x         y    z(sun)
  TRIAD sun-first  predicted       0.050     0.050     1.155
  TRIAD sun-first  simulated       0.050     0.050     1.149
  TRIAD mag-first  predicted       0.050     1.000     1.155
  TRIAD mag-first  simulated       0.051     1.001     1.149
  QUEST            predicted       0.050     0.050     1.155
  QUEST            simulated       0.050     0.050     1.149

Case 3: exact 180 deg rotation (QUEST singularity)
  sequential rotation about axis Some("y"), attitude error 2.23e-16 rad

Case 4: cost per solve, 6 observations (200000 solves, checksum 70306.747)
  QUEST    0.13 us
  q-method 1.01 us  (Jacobi eigen-solve)
```

Timings vary by machine.

## Takeaways

- Star tracker roll about the boresight is ~11x worse than cross-boresight,
  because roll is only seen through the few-degree spread of the stars. This
  is why trackers are often mounted in pairs at an angle.
- With one accurate and one coarse vector, TRIAD with the accurate vector
  first is essentially optimal (rotation about the Sun line is
  `sigma_mag / sin 60 deg` = 1.155 deg for both). Putting the magnetometer
  first throws away the sun sensor on one axis: 1.0 deg instead of 0.05 deg.
- QUEST gives the q-method answer about 8x faster here. Its only weak spot
  is 180 deg rotations, which the frame-selection step removes.
