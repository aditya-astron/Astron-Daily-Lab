# 025 - Polyphase DFT Channelizer for a Multi-Satellite Ground Receiver

**Domain:** Signal processing

## Problem

A ground-station SDR digitises an 8 Msps complex band that carries up to 16
CubeSat downlinks on a 500 kHz raster. Every link needs its own 500 ksps
baseband stream. Mixing, filtering and decimating each channel separately
costs about 14k real multiplies per input sample; split the band with a
single critically sampled polyphase DFT filter bank instead, keep 80 dB of
adjacent-channel rejection, and show that a weak link sitting between two
neighbours 30 dB hotter still demodulates cleanly.

## Approach

Channel k is centred at `f_k = k fs/M` (k > M/2 are negative frequencies).
The reference definition is mix, low-pass, decimate:

    y_k[m] = sum_l h[l] x[mM - l] exp(+j 2 pi k l / M)

With `l = rM + p`, this factors into M short branch filters at the low rate
followed by an inverse DFT across branches:

    h_p[r] = h[rM + p],   x_p[m] = x[mM - p]
    u_p[m] = sum_r h_p[r] x_p[m - r]
    y_k[m] = sum_p u_p[m] exp(+j 2 pi k p / M) = M * IFFT_p(u)[k]

**Prototype.** Kaiser-window low-pass, cutoff at half the spacing (250 kHz),
passband edge 200 kHz, stopband edge 300 kHz, 80 dB. The hand Kaiser estimate
`N = (A - 7.95)/(2.285 dw) + 1` gives 402 taps; the design rounds up to the
smallest odd length whose group delay `(N-1)/2` is a whole number of output
samples (417 taps, 208 = 13 x 16), then zero-pads to 432 = 27 x 16.

**Alias rule.** Critical sampling (output rate = spacing) folds anything
beyond +/-250 kHz back into the band. Energy past 300 kHz is 80 dB down before
it folds; energy in the 250-300 kHz transition is not, and lands in the outer
50 kHz of the neighbour. So each link must occupy at most +/-200 kHz. A
250 ksym/s QPSK link with RRC roll-off 0.35 occupies +/-169 kHz.

**Cost.** Per complex input sample: direct = `M (4 + 2N)` real multiplies,
polyphase = `(2N + 4 (M/2) log2 M) / M` = 62, a factor of 224.

**Scenario.** Five RRC-QPSK links (channels 2, 5, 11 at 0/0/-6 dB, channels
4 and 6 at +30 dB), matched RRC filter on each channel output, EVM after a
least-squares gain fit. No noise is added, so the EVM is the distortion floor
of the whole chain. The transmit RRC is 48 symbols long; at 12 symbols its
own truncation sidelobes (-53 dB) would dominate the result.

**Limits.** Critically sampled only (a 2x oversampled bank would relax the
transition band), ideal channel-centred carriers, no frequency offset or
timing recovery.

## Verification

`test_channelizer.py` checks:

- polyphase output equals the per-channel mix/FIR/decimate reference on
  random input (all 16 channels, relative error < 1e-10)
- prototype stopband >= 80 dB past 300 kHz, passband ripple within the Kaiser
  bound `20 log10((1+d)/(1-d))`, tap count against the hand estimate
- impulse response peaks exactly 13 output samples late with value h[208]
- tones at centre, +/-195 kHz and in negative-frequency channels: unity gain
  matching `|H(f)|` to 0.001 dB, every other channel below -80 dB
- alias rule: a tone at +280 kHz folds to -220 kHz at `|H(280 kHz)|`
  (-29 dB); one at +320 kHz folds below -80 dB
- multiply counts against the hand numbers (13888 and 62)
- QPSK scenario: all links EVM < -55 dB, the weak link between two +30 dB
  neighbours within 10 dB of an isolated link, empty channels < -80 dBc

## Run

    pip install numpy scipy
    python3 test_channelizer.py      # or: pytest
    python3 channelizer.py

## Sample output

```
fs = 8 Msps, M = 16, spacing 500 kHz, output 500 ksps per channel
prototype: 432 taps (27 per branch), group delay 13 output samples
  passband ripple 0.00165 dB (0-200 kHz), stopband 80.1 dB (>= 300 kHz)
real multiplies per input sample: direct 13888, polyphase 62.0 (224x fewer)

QPSK 250 ksym/s, RRC 0.35, 3000 symbols per link, no noise
 ch   centre kHz   level dB   EVM dB
  2       +1000      +0.0     -64.1
  5       +2500      +0.0     -57.3
 11       -2500      -6.0     -64.2
  4       +2000     +30.0     -64.8
  6       +3000     +30.0     -64.8
worst empty-channel power: -85.6 dB relative to channel 4
```
