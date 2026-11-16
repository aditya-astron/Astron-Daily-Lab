# 041 - CCSDS 121.0 Adaptive Rice Compression of Spacecraft Telemetry

**Domain:** Telemetry data compression

## Problem

A 3U CubeSat produces about 5.5 MB of sampled data a day: housekeeping
thermistors, a magnetometer and a 12-bit pushbroom imager. Its downlink
passes cannot carry all of that raw. CCSDS 121.0-B, the lossless "Rice"
coder flown on many missions, is cheap enough for a flight processor.
This project implements the full coder and decoder in C, checks that it
produces valid streams against libaec (the reference implementation), and
compares the coded rate with the entropy limit.

## Approach

**Preprocessor.** A unit-delay predictor `x_hat = x[i-1]`. The residual
`D = x - x_hat` is mapped to a non-negative integer using
`theta = min(x_hat, 2^n - 1 - x_hat)`:

```
delta = 2D            0 <= D <= theta
        2|D| - 1      -theta <= D < 0
        theta + |D|   otherwise (one-sided tail)
```

For any prediction this is a bijection onto `[0, 2^n - 1]`, so the mapped
residual never needs more than n bits.

**Entropy coder.** Each block of J residuals gets the shortest of:

| option | ID (n = 9..16) | cost |
|---|---|---|
| zero-block | `0000` + `0` | one FS code for the run length (1-4, ROS, >=5) |
| second extension | `0000` + `1` | FS(gamma) per pair, `gamma = (a+b)(a+b+1)/2 + b` |
| split-sample k | `k+1` | FS(delta >> k) for each sample, then each k-bit LSB field |
| no compression | `1111` | n bits per sample |

FS(m) is m zeros followed by a one. The coder tries every k (0 to 13 for
n = 16), so the option it picks is always the cheapest. Runs of zero
blocks never cross a 64-block segment. A run that reaches the end of its
segment and is 5 or more blocks long is coded as ROS ("remainder of
segment"). Every `rsi` blocks the first sample is sent raw as a
reference, right after the option ID, so each RSI can be decoded on its
own. The ID width is 3, 4 or 5 bits for n <= 8, 9..16 or > 16.

**Optimal k for geometric residuals.** If `P(delta = j) = (1-t) t^j`,
split-sample k costs `k + 1 + q/(1-q)` bits per sample with
`q = t^(2^k)`. Kiely (IPN PR 42-159, 2004) gives the minimiser in closed
form: `k* = max(0, 1 + floor(log2(ln(phi-1) / ln t)))`.

## Verification (`make test`)

`test_rice` runs 12 checks:

- The mapper is a bijection with an exact inverse for every (x, x_hat),
  checked exhaustively for n = 1..12. Spot values are checked for both the
  interleaved branch and the one-sided tail.
- 1280 random round trips cover n = 1..16, J = 8..64, RSI = 1..200, partial
  final RSIs, and rail-hugging and plateau data. Forced k = 0..13 round
  trips too.
- 16-bit white noise: 99.7 % of blocks go raw and the size never exceeds
  `blocks * (L + nJ)`.
- Constant data costs exactly 36 bits per 128-block RSI: two ROS codewords,
  one of them carrying the reference sample. That gives 90 bytes, the same
  as libaec.
- Zero-block runs of 1 to 12 blocks in the middle of a segment are coded
  without ROS and decode correctly.
- Second extension carries more than half the blocks on sparse +-1 data.
  The adaptive coder is never longer than the best fixed k.
- With k forced, the coded length matches `k + 1 + q/(1-q)` within 5 sigma
  in all 32 (theta, k) cases.
- Kiely's closed-form k* equals the brute-force argmin over 438 values of
  theta.
- A truncated stream returns an error at every cut point and never reads
  past the end of the buffer.

`test_libaec_interop.py` runs 112 cases: n = 8, 12, 14 and 16, four
(J, RSI) pairs, and signals that exercise every option, including
partial RSIs. In each case our stream is decoded by libaec, and libaec's
stream is decoded by ours. Both must reproduce the input exactly. Our
output is never longer than libaec's. It is the same length in every case and
byte-identical in 36.

## Run

```
make            # rice_demo, rice_tool, test_rice
make test       # C checks + libaec interop (pip install imagecodecs numpy)
make run
```

## Sample output

```
product                     n  samples    raw kB  coded kB    CR  b/smp     H0   %zb  %se %spl %raw
thermistor 1 Hz            12    86400     129.6      21.2  6.10  1.967  1.593   0.0 18.2 81.8  0.0  ok
magnetometer 10 Hz         16   864000    1728.0     493.6  3.50  4.571  4.239   0.0  0.0 100.0  0.0  ok
imager 2048 x 1024         12  2097152    3145.7    1770.5  1.78  6.754  6.464   0.0  0.0 100.0  0.0  ok
encrypted payload          16   262144     524.3     532.5  0.98 16.250 15.808   0.0  0.0  0.3 99.7  ok
                           worst-case expansion bound 1 + L/(nJ) = 1.0156

day total: 5527.6 kB raw -> 2817.8 kB coded (1.96:1)

White Gaussian noise about mid-scale (n = 16, J = 16): coded rate vs. entropy
   sigma    b/smp       H0     h(D)  rate-H0
       4    4.880    4.549    4.547    0.331
      64    8.864    8.541    8.547    0.322
    1024   12.859   12.484   12.547    0.376

Block size trade at sigma = 8:   J = 8: 6.074   16: 5.867   32: 5.769   64: 5.723 b/smp

Geometric residuals: theta 0.95 -> k* = 4, L = 5.786 vs H = 5.728 bits
```

How to read the results:

- **Gap to entropy.** The coded rate sits about 0.32 bit/sample above the
  residual entropy H0. Of that, 0.25 bit is the 4-bit ID spread over
  J = 16 samples. The rest is Rice-code redundancy.
- **Theory check.** For white noise with sigma >= 4, the measured H0
  matches the Gaussian differential entropy `0.5 log2(4 pi e sigma^2)`
  of the unit-delay residual.
- **Block size.** Larger blocks amortise the ID, but they adapt less well
  to data that changes along the stream.
- **Incompressible data.** Data that is already compressed or encrypted
  grows by exactly the ID overhead, 1/64 for this configuration.

## Limits

- Only the unit-delay predictor and unsigned samples, n <= 16. Signed data,
  external predictors and the restricted ID set for low n are not built.
- The 121.0-B-3 low-entropy extensions are not built either.
- The 2D predictor of CCSDS 123.0 would do better on the imager, because
  a pixel's neighbour in the previous row predicts it better than a 1D
  delay. The imager is included to show what a 1D coder leaves on the
  table.
