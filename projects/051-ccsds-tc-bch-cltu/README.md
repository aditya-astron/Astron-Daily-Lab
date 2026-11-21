# 051 - CCSDS Telecommand BCH(63,56) Codec and CLTU Link

**Domain:** Telemetry framing / Error correction

## Problem

Every command to a CCSDS spacecraft goes up as a CLTU: a 16-bit start
sequence `EB90`, a run of 8-octet BCH(63,56) codeblocks, and the tail
`C5C5C5C5C5C5C579` (CCSDS 231.0-B). The onboard decoder runs in one of two
modes: triple-error detection (accept only a zero syndrome) or single-error
correction. The questions here are the ones a TC link designer has to answer
with numbers, not folklore: how often each mode delivers, rejects or silently
miscorrects a block at a given BER; why the tail sequence ends a CLTU; and
what the probability of getting a whole frame through actually is.

## Approach

**The code.** `g(x) = x^7 + x^6 + x^2 + 1 = (x + 1)(x^6 + x + 1)`. With
`x^6 + x + 1` primitive, this is the (63,57) Hamming code cut down to its
even-weight half: `d = 4`, so SEC also detects every double error. Parity
bits are complemented and a zero filler bit pads the block to 64 bits.

Three encoders are implemented and must agree: the bit-serial LFSR of the
standard, a byte-wise table encoder (the 7-bit register left-aligned in a
byte is just CRC-8 with polynomial `0x8A`), and generic GF(2) division.

**Syndromes.** Column `i` of H is `x^i mod g`. Single errors give the 63
odd syndromes with nonzero `x^6+x+1` residue; every double error gives an
even nonzero syndrome. That leaves exactly one odd syndrome, `p(x) = 0x43`,
that no 1- or 2-bit error produces.

**Weight distribution, three ways.** MacWilliams from the 128-word dual
(simplex + complements: weights 0, 31, 32, 63), the closed form for the even
part of the Hamming enumerator
`A(z) = [(1+z)^n + n(1-z)(1-z^2)^((n-1)/2)]/(n+1)`, and `A_4` counted
directly from H columns (pairs with equal XOR). All in exact `i128`.

**Coset enumerator.** For any received word `t`,

    B_w = (1/|C_perp|) sum_{u in C_perp} (-1)^<u,t> K_w(wt u)

with Krawtchouk `K_w`. Because `<u_a, t> = <a, syndrome(t)>`, the coset
weights depend only on the syndrome. This gives the exact distance of the
tail from the code and the chance that noise turns it into a "valid" block.

**Exact BSC outcomes.** Radius-1 spheres around words at distance >= 4 are
disjoint, so the probability of landing in one is a plain sum over weights:

    P(sphere) = sum_w A_w [ P(w) + w P(w-1) + (n-w) P(w+1) ],  P(e) = p^e q^(n-e)

Detection mode uses radius 0. Nonzero codewords give the miscorrection
(undetected) probability.

**CLTU link.** TC frame (5-octet header, data, CRC-16-CCITT FECF),
randomized with `x^8+x^6+x^4+x^3+x^2+x+1`, split into codeblocks with `0x55`
fill, framed by start and tail, surrounded by idle `0x55`. The receiver
searches bit by bit for `EB90` (0 or 1 error allowed), decodes blocks until
one is rejected, derandomizes and checks length and FECF. Predicted delivery:

    P = P_start * P_c^(N-1) * P_last

where `P_last` is the probability the last block decodes to anything whose
frame bits are right: only 8 of its 56 info bits are frame here, the rest is
fill. That is the sphere probability over the *shortened* (55,48) code, whose
dual is the punctured TC dual, so MacWilliams handles it too. Leaving this
out under-predicts delivery by 4 sigma at BER 2e-2.

## Verification (`cargo test --release`)

- `g = (x+1)(x^6+x+1)`, `x` has order 63 mod `x^6+x+1`.
- LFSR, table and polynomial encoders agree on 20 000 blocks; all-zero data
  gives `00..00 FE`.
- `A_w`: MacWilliams = closed form, `A_4 = 9765 = n(n-1)(n-3)/24` also from
  column pairs, `sum A_w = 2^56`, dual weights as above.
- MacWilliams, the coset formula and the shortened distribution match brute
  enumeration on the (15,10) sibling code, as do all exact outcome
  probabilities (all 2^15 error patterns through a syndrome decoder).
- Decoder exhaustive to weight 3: all 63 singles corrected, all 1953 doubles
  rejected, all 39 711 triples rejected in detect mode; in SEC mode exactly
  `4 A_4 = 39 060` are miscorrected and the other 651 are the weight-3 words
  of syndrome `p(x)`.
- Tail: syndrome is `p(x)`, coset leader weight 3, rejected in both modes,
  still rejected with any single error, accepted (SEC) with every one of the
  1953 double errors.
- CRC-16 check value `0x29B1`; randomizer starts `FF 39 9E 5A 68` (CCSDS
  231.0-B), period 255, 128 ones.
- CLTU round trip for frames of 7 to 1024 octets; one flip in every block is
  fully repaired in SEC and kills the frame in detect mode.
- Monte Carlo decoder within 4 sigma of exact; small-BER leading terms
  `A_4 p^4`, `C(63,2) p^2`, `4 A_4 p^3`.

## How to run

    cargo run --release      # report (about 9 s)
    cargo test --release     # verification

Only the Rust standard library is used.

## Sample output

```
CCSDS TC code: g(x) = x^7 + x^6 + x^2 + 1 = (x + 1)(x^6 + x + 1), (n, k) = (63, 56)
d_min = 4, codewords = 2^56
weight distribution (MacWilliams from the 128-word dual, = closed form):
   A_4  = 9765
   A_6  = 1057224
   A_8  = 60544953
   A_30 = 13449656041565856
   A_32 = 14317376396958243
   dual weights: 0 x1, 31 x63, 32 x63, 63 x1

1) Codeblock outcomes on a BSC (exact, from the weight distribution)
    BER   |  detect-only: reject      undetected  |  SEC: correct+fixed  reject      miscorrect
      1e-6 |     6.300e-5       9.764e-21   |   0.999999998        1.953e-9   3.906e-14
      1e-5 |     6.298e-4       9.759e-17   |   0.999999805        1.952e-7   3.904e-11
      1e-4 |     6.281e-3       9.708e-13   |   0.999980549        1.941e-5    3.883e-8
      1e-3 |     6.109e-2        9.206e-9   |   0.998124663        1.839e-3    3.680e-5
      1e-2 |     4.690e-1        5.457e-5   |   0.868754525        1.094e-1    2.182e-2
      3e-2 |     8.518e-1        1.455e-3   |   0.432726565        3.652e-1    2.020e-1

2) Monte Carlo decoder vs exact (2e6 random codeblocks each)
   mode     BER    | outcome      exact        simulated    z
   Detect   0.03   | correct     1.467639e-1  1.465195e-1  -0.98
   Detect   0.03   | reject      8.517814e-1  8.519950e-1  +0.85
   Detect   0.03   | undetected  1.454735e-3  1.485500e-3  +1.14
   Correct  0.01   | correct     8.687545e-1  8.689205e-1  +0.70
   Correct  0.01   | reject      1.094296e-1  1.094405e-1  +0.05
   Correct  0.01   | undetected  2.181591e-2  2.163900e-2  -1.71
   Correct  0.03   | correct     4.327266e-1  4.325750e-1  -0.43
   Correct  0.03   | reject      3.652376e-1  3.653580e-1  +0.35
   Correct  0.03   | undetected  2.020359e-1  2.020670e-1  +0.11

3) Tail sequence C5C5C5C5C5C5C579
   syndrome = 0x43, filler bit = 1, coset leader weight = 3
   coset weights: B_3 = 651, B_4 = 0
    BER   |  P(tail accepted, SEC)   P(tail accepted, detect)
      1e-5 |       1.952e-7              6.506e-13
      1e-3 |       1.839e-3               6.132e-7
      1e-2 |       1.094e-1               3.624e-4

4) CLTU link: 71-octet TC frame -> 11 codeblocks, 8 octets of 0x55 idle either side
   start sequence vs every idle/start straddling window: min distance 5
   last block carries 8 frame bits + 48 fill bits
   config            BER   | delivered (exact)  delivered (sim)    z   | CRC-caught  undetected  false locks
   Correct start<=0 1e-3  |    0.96403          0.96391       -0.28  |       65           0          18
   Correct start<=0 5e-3  |    0.59079          0.59123       +0.40  |     4704           0          99
   Correct start<=1 5e-3  |    0.63829          0.63727       -0.95  |     5138           0         533
   Detect  start<=0 5e-3  |    0.02861          0.02830       -0.84  |        0           0          12
   Correct start<=1 2e-2  |    0.00773          0.00783       +0.55  |     5870           0        4435
   (200000 CLTUs per row)
```

Reading it: SEC buys about two orders of magnitude in block rejection at
BER 1e-4 (1.9e-5 vs 6.3e-3) and pays with miscorrections that grow as
`4 A_4 p^3` instead of `A_4 p^4`; at 1e-3 that is 3.7e-5 per block, which is
why the frame CRC is not optional in SEC mode (65 to 5870 frames per run were
caught by it, none slipped through). The tail is not merely "invalid": it
sits on the single syndrome that is three bit flips from the code, so no
single error can fake a block, but any double error makes it correctable.
That is harmless, because the next block is idle and the frame length field
discards anything extra. Allowing one error in the start sequence adds 4.7
points of delivery at 5e-3 and is safe against the idle pattern (nearest
straddling window is 5 bits away); the extra false locks it causes come from
searching random data after a missed start or a rejected first block, and none
of them produced a wrong frame.
