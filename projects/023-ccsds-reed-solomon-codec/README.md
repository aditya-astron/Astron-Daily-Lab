# 023 - CCSDS Reed-Solomon (255,223) Codec with Errors-and-Erasures Decoding

**Domain:** Error correction / telemetry

## Problem

Project 002 framed CCSDS telemetry with a CRC, which can only reject a damaged
frame. The outer code that actually repairs frames on CCSDS links is the
Reed-Solomon (255,223) code: 32 parity bytes per 223 data bytes, correcting any
16 wrong bytes in a codeword. This project implements the full codec as the
standard specifies it (field, generator, dual-basis symbols, virtual fill,
interleaving) plus an errors-and-erasures decoder, and checks it from three
directions: algebra, an independent library, and channel statistics.

## Approach

Code definition (CCSDS 131.0-B):

```
GF(2^8):    F(x) = x^8 + x^7 + x^2 + x + 1,  alpha a root
generator:  g(x) = prod_{j=112}^{143} (x - beta^j),  beta = alpha^11
symbols:    Berlekamp dual basis on the wire
```

The roots pair up as beta^(112+j) * beta^(143-j) = 1, so g(x) is palindromic,
which is why CCSDS picked these exponents (half the multipliers in hardware).

- **Encoder:** systematic LFSR division, parity = m(x) x^32 mod g(x).
- **Syndromes:** S_j = r(beta^(112+j)), j = 0..31.
- **Erasures:** the erasure locator Gamma(x) = prod(1 + X_k x) seeds
  Berlekamp-Massey, which then only has to find the unknown errors. Any
  combination with 2e + s <= 32 is decodable.
- **Chien search:** roots of Lambda(x) give the bad positions. Because the code
  is full length, every nonzero field element is a position, so finding fewer
  roots than deg Lambda means the word is out of range: a decoder failure
  rather than a miscorrection.
- **Forney:** for first root beta^112,
  `Y = X^(1-112) Omega(X^-1) / Lambda'(X^-1)`, Omega = S Lambda mod x^32.
- **Dual basis:** the CCSDS T-matrix maps conventional bytes to the wire. The
  test proves bit (7-m) of the wire symbol is Tr(x alpha^(117 m)), i.e. the
  matrix is the trace dual of the basis {alpha^(117 m)}.
- **Shortening:** virtual zero fill in front of the data, as in CCSDS.
- **Interleaving:** depth I puts byte m in codeword m mod I, so any channel
  burst of 16 I bytes is corrected.

Channel analysis in `rs_link.c`: with independent symbol errors at rate ps
(what an ideal interleaver delivers) a bounded-distance decoder fails with

```
P_W = sum_{i=17}^{255} C(255,i) ps^i (1-ps)^(255-i)
```

and the decoded bit error rate follows by weighting each failure by i/255
wrong symbols and 128/255 wrong bits per wrong symbol.

## Files

| File | Purpose |
|------|---------|
| `rs.h`, `rs.c` | codec: field tables, encoder, errors-and-erasures decoder, dual basis, shortening, interleaving |
| `test_rs.c` | 470 checks: field and generator algebra, dual basis, linearity, every e <= 16, erasure mixes up to s = 32, refusal past capacity, shortened and interleaved frames |
| `verify_galois.py` | cross-check of generator, encoder and decoder against the `galois` library via ctypes |
| `rs_link.c` | Monte Carlo vs. bounded-distance formula, coding gain on BPSK, burst tolerance vs. depth |

## How to run

```
make test            # unit checks
make run             # link analysis
pip install galois   # optional, then:
make verify          # independent cross-check
```

## Sample output

```
$ ./test_rs
470/470 checks passed

$ ./rs_link
CCSDS RS(255,223), t = 16, rate = 0.8745
g(x) coefficients, x^32 first (conventional basis):
  01 5b 7f 56 10 1e 0d eb 61 a5 08 2a 36 56 ab 20 71 20 ab 56 36 2a 08 a5 61 eb 0d 1e 10 56 7f 5b 01

Independent symbol errors, 20000 words per point
   ps     P_W formula   P_W sim     fails  miscorr   z
  0.030  2.0161e-03    2.2500e-03      45       0   +0.74
  0.035  9.0593e-03    7.9500e-03     159       0   -1.66
  0.040  2.8818e-02    2.8900e-02     578       0   +0.07
  0.045  7.0624e-02    6.9000e-02    1380       0   -0.90
  0.050  1.4139e-01    1.4240e-01    2848       0   +0.41

Hard-decision BPSK, ideal symbol interleaving, target BER 1e-5
  uncoded Eb/N0      9.59 dB
  RS(255,223) Eb/N0  6.29 dB   (channel ps = 2.509e-02)
  coding gain        3.29 dB

Interleaving depth vs. longest correctable channel burst
  I   frame bytes  burst bytes  burst bits (guaranteed)
  1     255          16           121
  ...
  5    1275          80           633
  8    2040         128          1017

simulation agrees with theory

$ make verify
galois cross-check: OK
```

## Notes

- Every simulated failure past 16 errors was detected; none miscorrected.
  For RS(255,223) the miscorrection probability is roughly 1/t! = 1/16!,
  which is why the CRC behind it in project 002 almost never fires.
- The 3.3 dB here is the RS code alone on hard decisions. On real CCSDS links
  it sits outside the K=7 convolutional code of project 008, whose bursty
  Viterbi output is exactly why the interleaver exists.
