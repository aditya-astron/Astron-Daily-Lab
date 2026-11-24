# 062 - Sinc^N Decimation of a Second-Order Sigma-Delta Current-Sense Stream

**Domain:** Signal processing / Embedded firmware

## Problem

A PCDU shunt channel uses an isolated second-order 1-bit sigma-delta
modulator clocked at 10 MHz. The bitstream is decimated in firmware (or the
MCU's filter block) by a sinc^N / CIC filter. The current loop wants the
highest output rate and lowest latency that still gives 80 dB SNR on a
half-scale sine. Which sinc order N and decimation ratio R, how wide must the
registers be, and how many bits can Hogenauer pruning drop?

C11 for the modulator and the fixed-point decimator, Python (numpy, ctypes)
for the study.

## Approach

**Modulator** (`cic.c`, `sdm2_run`): NTF = (1 - z^-1)^2, STF = z^-2.

```
y[n] = sign(s2[n]);  s1[n+1] = s1[n] + x[n] - y[n];  s2[n+1] = s2[n] + s1[n] - 2 y[n]
```

With e = y - s2 this is exactly Y = z^-2 X + (1 - z^-1)^2 E.

**Decimator** (`cic.c`): N integrators, decimate by R, N combs of delay M, all
in unsigned registers that wrap modulo 2^width. Wraparound is harmless when
the width covers the true output range:

```
Bmax = ceil(N log2(RM)) + Bin
```

Pruning discards B_j LSBs at the input of stage j. A logical shift of the
unsigned residue is floor division modulo the narrower width, so pruning
commutes with wraparound. Hogenauer's rule:

```
F_j^2 = sum_k h_j[k]^2,  H_j = (1 - z^-RM)^N / (1 - z^-1)^(N-j+1)  (j <= N)
                         H_j = (1 - z^-M)^(2N+1-j)                  (combs)
B_j   = floor(-log2 F_j + log2 sigma_T + 0.5 log2(6/N)),  sigma_T^2 = 2^(2 B_2N+1) / 12
```

The error model uses exact discrete truncation moments, mean
-(2^Bj - 2^Bj-1)/2 and variance (4^Bj - 4^Bj-1)/12; only stage 1 and the
output have nonzero DC gain, so only they bias the output.

**Noise vs. R** (`decimation_study.py`). If e is white, the output noise is
sigma_e^2 sum g^2 / (RM)^2N with G = (1 - z^-1)^2 (boxcar_RM)^N:

| N | G(z) | sum g^2 | SNR growth |
|---|------|---------|-----------|
| 1 | (1 - z^-1)(1 - z^-R) | 4 | 6 dB/oct |
| 2 | (1 - z^-R)^2 | 6 | 12 dB/oct |
| 3 | (1 - z^-R)^2 boxcar_R | 6R | 15 dB/oct |

For N <= L the shaped noise near fs/2 folds into the output band; N = L + 1
is the first order that gets the full 15 dB/octave.

## Results

```
1. SNR vs decimation ratio: simulation / linear-model prediction [dB]
        R               N=1               N=2               N=3               N=4
       16     15.7 / 11.7       34.9 / 34.0       45.0 / 46.0       48.5 / 49.5
       64     27.4 / 23.7       58.7 / 58.1       75.1 / 76.1       78.7 / 79.6
      256     39.6 / 35.7       82.8 / 82.1      104.8 / 106.2     108.2 / 109.7
   slope dB/octave:                  6.0              12.0              14.9              14.9
   quantizer error variance sigma_e^2 = 0.545 (uniform-error model 1/3 for +/-1 levels)

2. Smallest R for SNR >= 80 dB (half-scale sine)
     N  R model  R sim     out rate    pred dB   sim dB      f-3dB   delay [us]   Bmax
     2      227    227    44.1 kS/s       80.0     80.7    14.0 kHz      22.6     18
     3       77     81   123.5 kS/s       81.2     80.4    32.3 kHz      12.0     22
     4       66     69   144.9 kS/s       81.3     80.3    32.9 kHz      13.6     27

3. Hogenauer pruning at N = 3, R = 81 (full width 22 bits)
    bout discarded B_1..B_7              prune noise  pred loss   sim loss
      18 0 0 0 0 1 2 4                     -100.8 dB       0.30       0.30
      17 0 0 0 1 2 3 5                      -94.6 dB       1.15       1.19
      16 0 0 0 2 3 4 6                      -88.5 dB       3.49       3.34

Design: sinc^3, R = 81 -> 123.5 kS/s, 12.0 us group delay, 18-bit output word,
        registers 22/22/22/22/21/20 bits, simulated SNR 80.1 dB
```

Takeaways:

- sinc^3 gives 2.8x the output rate of sinc^2 at half the latency for the same
  SNR. sinc^4 buys another 17 % rate but costs a stage, 5 more register bits
  and more group delay than sinc^3 at its own ratio.
- The 1-bit quantizer error is not uniform (sigma_e^2 = 0.545, not 1/3) and
  not quite white, so the white model is ~1 dB optimistic for N >= 3; R is
  stepped up from 77 to 81 by simulation.
- Pruning to an 18-bit word costs 0.3 dB; 16 bits would cost 3.3 dB and miss
  the requirement. The first three integrators cannot be pruned at all.
- A +/-1 stream never reaches -2, so 21 bits would do for this input; the 22
  bit formula covers any 2-bit input.

## Verification

`test_cic.c` (19 checks):

| Check | Result |
|---|---|
| Bit-exact vs a 128-bit direct-form boxcar cascade, 4 configs incl. M = 2, 16-bit input | exact |
| Bmax tight: constant -2 correct at 22 bits, wraps at 21; +/-1 extremes fit 21, wrap at 20 | as predicted |
| Hogenauer F_j^2 closed form vs brute-force impulse responses, 3 configs | exact |
| Pruned N=4, R=32 to 16 bits: error variance and mean vs model | 0.08 %, -23.99 vs -24.00 |
| Modulator residual of Y = z^-2 X + (1 - z^-1)^2 E | 1e-15 |
| Bitstream mean for DC 0.3 | 0.299999 |

`test_decimation.py` (13 checks): C decimator equals a numpy boxcar cascade on
a real bitstream; output splits exactly into CIC(signal) + CIC(NTF e);
sum g^2 closed forms 4, 6, 6R; measured slopes within 0.5 dB/oct of 6, 12,
15; white model within 1.5 dB for N = 2..4; sinc^3 droop equals (2/pi)^3 at
output Nyquist and f-3dB matches the sinc^3 root; N = 3, R = 81 meets 80 dB
while N = 2 at the same R misses by over 15 dB; pruned 18-bit design meets
80 dB with loss matching the Hogenauer prediction to 0.1 dB.

## Run

```
make test     # C checks, then pytest
make run      # design study
```
