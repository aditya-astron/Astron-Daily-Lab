# 047 - GNSS snapshot positioning with RAIM fault detection and exclusion

A single-frequency GNSS receiver has to tell the user when its own fix can't be
trusted, with no ground monitoring: one satellite with a bad clock or
ephemeris can pull the solution hundreds of metres. Receiver autonomous
integrity monitoring (RAIM) uses the redundancy beyond four satellites to
detect that, bound the horizontal error it could hide (HPL), and drop the
faulty satellite. This project builds the fix, the residual test, the
protection levels and the exclusion step, checks every statistical claim by
Monte Carlo, and runs a 24 h availability study for a receiver in Mumbai on a
24-satellite GPS-like constellation.

## Approach

**Fix.** Pseudoranges `rho_i = |R3(we tau_i) r_s,i - r| + b + e_i` with the
Sagnac rotation (light time iterated). Gauss-Newton on `(r, b)` from Earth's
centre; ENU geometry matrix rows `[-cos el sin az, -cos el cos az, -sin el, 1]`.

**Residual test.** With `A = (H^T H)^-1 H^T` and `S = I - H A`:

```
fault free:      SSE / sigma^2 ~ chi2(N-4)              alarm if > T = chi2_inv(1-Pfa, N-4)
bias B on sat i: SSE / sigma^2 ~ ncchi2(N-4, lambda),   lambda = B^2 S_ii / sigma^2
horizontal shift = B |A_h,i| = slope_i sigma sqrt(lambda),  slope_i = |A_h,i| / sqrt(S_ii)
HPL = max_i slope_i * sigma * sqrt(lambda_det),  P(ncchi2(N-4, lambda_det) < T) = Pmd
```

The chi-square quantile works on whichever incomplete-gamma tail is small so
it stays accurate at `1 - 1e-7`; the noncentral CDF is the Poisson mixture of
central chi-squares summed outward from the Poisson mode.

**Exclusion.** On alarm (N >= 6), re-solve every N-1 subset and keep the one
that passes its own `N-5` dof test with the smallest statistic. For
availability, FDE counts as available when every N-1 subset has a detection
HPL under the HAL (the usual simplified criterion).

**Constellation.** Six planes at 55 deg, four unevenly spaced slots each,
circular orbits at a = 26559.7 km, mask 5 deg. Receiver at 19.076 N 72.878 E.

## Run

```
make run     # study
make test    # 43 C++ checks + scipy cross-check of the statistics (~75 s)
```

## Sample output (abridged)

```
1. Snapshot: 7 satellites, noiseless fix in 6 iterations, error 2e-10 m
   GDOP 2.25  PDOP 1.99  HDOP 1.09  VDOP 1.66
   sigma 5 m, Pfa 1e-5, Pmd 1e-3: T = 25.902 (3 dof), lambda_det = 64.381
   worst slope D3 1.094 -> HPL 43.9 m, VPL 61.9 m (sigma 33.3 m: HPL 292 m)

2. Bias ramp on D3
     bias m   SSE/s^2  horiz m  alarm  excluded horiz after m
         25      4.82    22.58     no         -        22.58
         50     27.96    24.94    yes        D3         3.73
        150    291.80    87.66    yes        D3         7.99
   72.7 m bias is missed 1e-3 of the time; its mean horizontal error is the HPL

3. 24 h, HAL 556 m (0.3 nmi NPA)
   sigma  5.0 m: HPL median 32 m, max 78 m
      all 24 up:        detection 100.00 %   FDE  98.96 %
      worst one out A1: detection  99.17 %   FDE  92.99 %
   sigma 33.3 m (SA-era DO-208 value): HPL median 215 m, max 518 m
      all 24 up:        detection 100.00 %   FDE  61.25 %
      worst one out F2: detection  86.39 %   FDE  49.51 %
```

With today's ranging error RAIM detection is essentially always there for
non-precision approach; exclusion is the part that thins out, and under the
old 33 m assumption a single missing satellite costs over three hours of
detection a day, which is why SA-era receivers needed baro-aiding and
prediction tools.

## Checks (`test_raim.cpp`, `verify_stats.py`)

- chi2 quantiles against table values and the 2-dof closed form; 1-dof
  noncentral CDF against its normal closed form; 128 quantiles and CDFs
  against scipy to 4e-15 relative and 1e-14 absolute.
- Noiseless fix from Earth's centre at 300 random sites and times to 1e-8 m
  within 6 iterations; ignoring Sagnac costs 33 m; solver ENU geometry matches
  the az/el construction to 6e-16.
- DOP of a zenith satellite plus a ring of n at elevation e against the closed
  form `HDOP^2 = 4/(n cos^2 e)`, `VDOP^2 = (n+1)/(n(1-sin e)^2)`,
  `TDOP^2 = (n sin^2 e + 1)/(n(1-sin e)^2)` for 24 geometries to 6e-14.
- trace(S) = N-4; a 100 m bias moves the nonlinear fix by `B A[:,i]` and gives
  `SSE = B^2 S_ii` to 1e-4.
- 100k fault-free epochs: ENU and clock variances match `sigma^2 (H^T H)^-1`
  within 2 %, SSE mean and variance match chi2(3), false-alarm rate matches Pfa.
- 40k epochs per satellite at `lambda_det`: missed-detection rate equals Pmd for
  every satellite, and the worst satellite's mean horizontal error equals HPL.
- Integrity over a 10-160 m bias sweep on all 7 satellites: the rate of
  "horizontal error above HPL and no alarm" never exceeds Pmd.
- 35k epochs with a 150 m fault: correct exclusion every time, no wrong ones.
