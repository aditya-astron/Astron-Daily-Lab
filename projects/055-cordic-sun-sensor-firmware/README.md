# 055 - Fixed-point CORDIC for a coarse sun sensor

**Domain:** Embedded firmware  **Language:** C11 (integer-only firmware, host-side test and study)

## Problem

A coarse sun sensor on a small spacecraft runs on an MCU with no FPU, and
safe-mode code should not pull in libm. The sensor has four cosine-law
photodiodes, and the firmware has to turn their ADC counts into a sun
direction. That takes an `atan2`, a vector magnitude and a second `atan2`.
The question is how many CORDIC iterations the firmware needs before its
arithmetic stops adding to the error the 12-bit ADC already has, and how to
prove the integer code stays within a known error bound for every input.

## Approach

**CORDIC.** Each micro-rotation is a shift and two adds:

    x' = x - d (y >> i),   y' = y + d (x >> i),   z' = z - d atan(2^-i)

Each step stretches the vector by `k_i = sqrt(1 + 2^-2i)`. After N steps the
gain is `K_N` (which tends to 1.6467602581), and one Q31 multiply removes it.
Angles are BAM (`2^32` = one turn), so wraparound costs nothing.

- *Rotation mode* (`cordic_sincos`) folds |angle| > 90 deg by pi and starts
  from `x = 1/K_N` in Q30, which gives cos and sin in Q2.30.
- *Vectoring mode* (`cordic_atan2_mag`) block-normalises the input so
  max(|x|,|y|) is in [2^28, 2^29). That leaves headroom for `sqrt2*K_N < 2.33`
  and keeps the error fixed relative to the vector length at any input scale.
  It then folds the left half-plane and drives y to zero. It handles any
  int32 input, including `INT32_MIN`.

**Error bounds** (in `common.h`, enforced by the test):

- Residual angle after N steps: `|z_N| <= atan(2^-(N-1))`
- Table rounding: `N/2` BAM
- Datapath: N-1 truncating shifts, each under 1 LSB per component and
  amplified by at most `K_N`, give `sqrt2 (N-1) K_N` LSB. In vectoring mode
  this is relative to `K_N r` with `r >= 2^28`.

**Sun sensor.** The cell normals are tilted by `t0 = 30 deg` from boresight
toward +x, -x, +y and -y. With `I_k = G n_k . s`:

    Dx = 2G sin t0 s_x,  Dy = 2G sin t0 s_y,  S = 4G cos t0 s_z
    azimuth = atan2(Dy, Dx),   off-axis = atan2(2 cot t0 |D|, S)

The gain G cancels, so cell ageing and solar distance drop out. The first
vectoring pass gives the azimuth and `|D|`. `|D|` is scaled by `2 cot t0`
(Q4.28) and fed to the second pass. Counts are lifted to Q12 so `|D|` keeps
fractional bits between the two passes. The firmware needs no square root,
no divide and no libm call.

## Verification

`test_cordic` checks:

- Both ROM tables are regenerated from libm and match bit for bit, and
  `K_inf = 1.6467602581210656`.
- Rotation mode: the worst error over 65536 angles is within the analytic
  bound for every N from 4 to 30. While the residual angle term dominates
  (N = 6..20), each extra iteration halves the error.
- Vectoring mode: angle and magnitude stay within the bounds over 200k random
  vectors spanning 1 to 2^31 in scale, plus the edge cases: the axes,
  `INT32_MIN` pairs, the zero vector, and 3-4-5 (which gives exactly 5).
- Sun sensor, hand case: the sun 40 deg off axis with G = 4000 gives counts
  {3939, 1368, 2654, 2654}, and `atan(2 sqrt3 * 2571 / 10615) = 39.9994 deg`.
  The solver matches this to 1e-5 deg. At boresight the solver returns 0, and
  `2 cot 30 deg` in Q28 is 929887697.
- Sun sensor, sweep: over 20k sun directions, the firmware agrees with a
  double-precision solve on the same counts.

It also runs clean under GCC `-fsanitize=undefined,address` and builds with
Clang. Running under the sanitizer caught a left shift of a negative count
difference, which is now a multiply.

**Hand check of the ADC floor.** Near boresight,
`d(off) = 2 cot t0 dD / S`, with `S = 4 * 4000 cos 30 = 13856`. For 1 count
of noise per cell, `sigma_D = sqrt2`, which gives 0.020 deg per axis and
0.028 deg RMS over both axes. The simulation gives 0.0274 deg. Quantisation
alone (`1/sqrt12` count) predicts 0.0081 deg; the simulation gives 0.0076 deg.

## Run

    make test     # builds and runs test_cordic
    make run      # accuracy table and sun-sensor error budget

## Sample output

```
CORDIC accuracy vs iterations (worst case over 65536 angles / 100000 vectors)
  N   sincos err   bound      atan2 err    bound      |r| rel err  bound
   8  7.812e-03  7.812e-03  7.812e-03  7.812e-03  3.051e-05  3.058e-05
  12  4.883e-04  4.883e-04  4.883e-04  4.885e-04  1.258e-07  2.036e-07
  16  3.052e-05  3.056e-05  3.052e-05  3.074e-05  1.381e-08  1.059e-07
  24  1.279e-07  1.874e-07  1.347e-07  4.303e-07  2.352e-08  1.475e-07
  30  1.637e-08  8.747e-08  2.503e-08  3.806e-07  3.523e-08  1.791e-07

Sun sensor pointing error, t0 = 30 deg, G = 4000 counts, 12-bit ADC, sigma = 1.0 count, FOV 55 deg
  solver            RMS [deg]   max [deg]   RMS over floor
  double precision  0.02741     0.10276     -
  CORDIC N=10       0.08024     0.19543    +192.80 %
  CORDIC N=12       0.03328     0.11134     +21.43 %
  CORDIC N=13       0.02898     0.10035      +5.74 %
  CORDIC N=14       0.02781     0.10220      +1.49 %
  CORDIC N=16       0.02743     0.10333      +0.09 %
```

## Findings

- Up to about N = 20, the worst-case error equals `atan(2^-(N-1))` and the
  datapath adds nothing measurable. With 32-bit words, the practical floor
  is about 2e-8 rad (sincos) and 3e-8 rad (atan2) at N = 28 to 30. The
  analytic datapath bound is conservative by roughly 5 to 15x there.
- For this sensor, N = 14 adds 1.5 % RMS to the 1-count noise floor and
  N = 16 adds 0.1 %. In a noiseless, quantisation-only case, N = 16 is
  needed to get within about 1 % of the floor. **N = 16** is the choice: two
  passes of 16 shift-add steps each, against 6 or 8 for a pass that would
  visibly degrade the sensor.
- Block normalisation matters more than the iteration count when inputs are
  small. Without it, 12-bit counts would use 12 of the 29 available bits,
  and the datapath floor would rise by about 2^17.
