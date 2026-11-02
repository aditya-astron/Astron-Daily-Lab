# 008 - CCSDS K=7 convolutional code with soft-decision Viterbi decoding

## Problem

The workhorse forward error correction on deep-space and LEO downlinks for
decades has been the CCSDS 131.0-B rate-1/2, constraint-length-7
convolutional code (generators 171 and 133 octal, second symbol inverted).
The goal is to build the encoder and a maximum-likelihood Viterbi decoder,
then show what the code actually buys on a BPSK/AWGN link: how much Eb/N0 it
saves, and how much of that is lost if the receiver only passes hard or 3-bit
decisions to the decoder.

## Approach

**Encoder.** A 6-bit shift register; each input bit produces two symbols
`c1 = parity(reg & 171)`, `c2 = parity(reg & 133) XOR 1`. Six zero tail bits
return the encoder to state 0 so the block can be decoded as a terminated
trellis.

**Viterbi decoder.** 64 states, add-compare-select per step, one survivor bit
per state packed into a 64-bit word. The branch metric is the correlation
`r1*e1 + r2*e2` with `e = +/-1`, which is the AWGN log-likelihood up to
constants. Feeding it `+/-1` turns it into a Hamming-distance decoder, and
feeding it 8 quantizer levels gives the classic 3-bit soft decoder, so one
routine covers all three receivers. Traceback starts at state 0 after the tail.

**Distance spectrum.** All paths that leave state 0 and first remerge are
enumerated by output weight, accumulating the number of events `a_d` and the
information-bit weight `B_d`. The code has no zero-weight cycle, so the
search terminates once every live path exceeds `d_max`.

**Performance bounds.**

- uncoded BPSK: `Pb = Q(sqrt(2 Eb/N0))`
- soft decision union bound: `Pb <= sum_d B_d Q(sqrt(2 d R Eb/N0))`
- hard decision (BSC, `p = Q(sqrt(2 R Eb/N0))`): `Pb <= sum_d B_d P_d`,
  with `P_d` the exact probability that more than `d/2` of `d` symbols flip
  (ties split evenly)
- asymptotic soft coding gain: `10 log10(R d_free) = 6.99 dB`

## Verification (`make test`)

- impulse response reproduces the 171/133 taps with the CCSDS inversion
- `d_free = 10` and `a_d`, `B_d` for d = 10..20 match the published table
  (11/36, 38/211, 193/1404, 1331/11633, 7275/77433, 40406/502690)
- noiseless 5000-bit block decodes exactly
- 20000 random blocks with 1 to 4 symbol errors (scattered and clustered) all
  correct, as guaranteed by `d_free = 10`; 6 errors placed on the minimum
  weight event do decode to the wrong codeword
- uncoded channel BER at 4 dB matches theory within 4 sigma
- coded soft BER at 4 dB sits below and within 3x of the union bound; soft
  beats 3-bit beats hard; hard is under its BSC bound

## How to run

```
make test     # unit and statistical checks, ~6 s
make run      # BER sweep, ~1 min; ./ber_sim N sets bits per point
```

## Sample output

```
d_free = 10, asymptotic soft coding gain = 10log10(R*dfree) = 6.99 dB

Eb/N0  uncoded    hard       3-bit      soft       union bound
 2.0   3.751e-02  1.117e-01  7.217e-03  5.405e-03  2.846e-02
 3.0   2.288e-02  3.411e-02  9.727e-04  3.886e-04  6.814e-04
 4.0   1.250e-02  5.000e-03  3.361e-05  1.737e-05  1.869e-05
 4.5   8.794e-03  1.650e-03  8.871e-06  1.624e-06  3.032e-06
 5.0   5.954e-03  5.648e-04  8.746e-07  8.746e-07  4.427e-07
```

Soft decoding reaches 1e-5 near 4.2 dB versus 9.6 dB uncoded, about 5.4 dB of
real gain, consistent with the usual 4.4 dB figure quoted for this code. The
3-bit quantizer costs roughly 0.2 to 0.3 dB, while hard decisions throw away
about 2 dB and are worse than no coding below about 3.5 dB. Points with only a
handful of error events (5.0 dB and up at 8e6 bits) carry large statistical
uncertainty; the union bound is the better guide there.
