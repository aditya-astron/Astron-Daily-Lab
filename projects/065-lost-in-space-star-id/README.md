# 065 - Lost-in-space star identification (k-vector + Pyramid)

A star tracker that wakes up tumbling, or after a reset, has no attitude prior.
It must name the stars in one frame from geometry alone, quickly, and must not
lock onto a wrong solution when the frame contains hot pixels, debris glints or
missing stars. This project builds that pipeline for a 12 x 12 deg, 1024 px
tracker and checks every stage against a closed form or a brute-force reference.

## Setup

- Catalog: 5000 stars drawn isotropically (roughly the count to visual
  magnitude 6). Synthetic, so the isotropic statistics below are exact
  references; a real catalog only changes the density map.
- Camera: pinhole, f = 512 / tan 6 deg = 4871.4 px, IFOV 42.3 arcsec/px,
  centroid noise 0.1 px (4.23 arcsec). Matching tolerance on a pair
  separation eps = 4 sqrt(2) sigma = 24 arcsec.

## Approach

**Pair catalog.** Every pair closer than the FOV diagonal (16.91 deg), sorted
by cos(separation). For an isotropic sky the count is
`N(N-1)/2 * (1 - cos d_max)/2`, and the separation pdf is `sin(t)/2`, so a
`t +- eps` query should return `N(N-1)/2 * (cos(t-eps) - cos(t+eps))/2` pairs.

**k-vector (Mortari 1997).** Lay a straight line `z(j) = q + m j` across the
sorted values and store `K[j] = #{y_i <= z(j)}`. A range `[ya, yb]` maps to
`K[floor((ya-q)/m) - 1] .. K[ceil((yb-q)/m)]` in O(1); the superset is about
one bin wider at each end and is trimmed by comparison.

**Pyramid (Mortari et al. 2004).** Pick observed triangles in the order
`dj, dk, i` (spacing outermost) so a spurious star is dropped quickly instead
of being retried in every triangle. For each side, query the k-vector; join
the three candidate lists through a neighbour map and a pair hash; keep only
catalog triangles with the same handedness (sign of `b_i . (b_j x b_k)`).
Near-collinear triangles are skipped because their handedness is within the
noise. A triangle is accepted only if a fourth observed star is matched by
exactly one catalog star at all three observed distances, and only one
candidate triangle is confirmed this way.

**Attitude.** Davenport q-method on the confirmed four, then every other
detection is a nearest-neighbour lookup within `2 eps`; a catalog star claimed
by two detections is dropped from both. Final q-method on all matches. The
expected error is Shuster's QUEST covariance `P = sigma^2 [sum(I - b b^T)]^-1`.
By hand, cross-boresight `sigma/sqrt(n) = 4.23/sqrt(17.4) = 1.0 arcsec`; roll
`sigma/(sqrt(n) rho_rms)` with `rho_rms = sqrt(2/3) * 6 deg = 0.0855 rad` gives
about 12 arcsec.

## Results (sample_output.txt)

| Check | Measured | Reference |
|---|---|---|
| Pairs within 16.91 deg | 271 016 | 270 149 (isotropic, +0.32 %) |
| Stars per frame | 17.34 | N Omega/4pi = 17.39, Omega = 4 asin(sin^2 6 deg) |
| k-vector vs binary search, 200k queries | 0 mismatches | identical index ranges |
| Candidates per side query, 5 bands | within 0.6 % | `sin(t)/2` pair density |
| Attitude RMS x/y/z at 0.1 px | 1.07 / 1.07 / 13.0 arcsec | QUEST 1.08 / 1.08 / 13.1 |
| Same, 0.05 to 1.0 px | within 1 % per axis | QUEST |

Robustness at 0.1 px, 2000 frames per row:

| Spikes | Drop-out | Correct | Wrong | No ID | ms/frame |
|---|---|---|---|---|---|
| 0 | 0 | 100 % | 0 | 0 | 0.2 |
| 10 | 0 | 100 % | 0 | 0 | 0.6 |
| 10 | 30 % | 99.6 % | 0 | 8 | 0.8 |

No frame in any run produced a wrong identification; the cost of spikes and
drop-outs shows up as a rare refusal, which is the failure mode you want.
The k-vector query takes about 60 ns against about 350 ns for two binary
searches (timings vary by machine). "Wrong" is judged by the labels, not an
attitude threshold: with few or clustered stars a correct solution can have a
roll error of an arcminute or more.

## Run

```
make test   # k-vector exactness, q-method, counts, identification and covariance checks
make run    # full study, about 12 s
```

Requires only a C++17 compiler.

## Tests

- k-vector against `lower_bound`/`upper_bound` on 100k queries with
  duplicates, clusters and exact endpoint hits.
- q-method: a 90 deg hand case and 200 random rotations to 1e-12.
- Pair list against brute force for 50 stars; FOV solid angle against
  focal-plane quadrature.
- 300 noiseless frames: every label right, attitude to 1e-9.
- 1500 frames with 0.1 px noise, 5 spikes, 20 % drop-out: zero wrong IDs,
  at least 99.5 % identified, per-axis RMS within 10 % of QUEST.
- 300 frames of spikes only: no solution returned.
