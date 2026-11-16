# 043 - Hsiao (72,64) SECDED EDAC and memory scrub interval design

## Problem

Flight memory sits behind an EDAC that corrects one flipped bit per 64-bit
word and detects two. Radiation keeps flipping bits, so a background scrubber
must read and rewrite every word often enough that a second upset rarely lands
in a word before the first one is cleaned. How often is "often enough", and
what does scrubbing not fix? This project builds the codec, proves its
properties exhaustively, derives the exact word-failure probability between
scrubs (including two-cell multi-cell upsets, MCU), checks that model against
a Monte Carlo of the real decoder, and sizes the scrub period for a 256 MiB
payload memory in GEO.

## Approach

**Code.** Hsiao odd-weight-column construction: the 64 data columns of H are
all 56 weight-3 columns of length 8 plus the 8 rotations of `00011111`;
check bit i has column `1<<i`. Every row then has 26 data ones (balanced XOR
trees). Odd columns mean a single error gives an odd syndrome and a double
error an even nonzero one, so SEC-DED needs no overall parity bit. Decoding
is a 256-entry syndrome table.

**Word failure between scrubs.** Each cell sees upset events at rate
lambda; a fraction f also flip the next physical cell. A bit's state after T
depends only on the parity of each source's Poisson count,
P(odd) = (1 - e^{-2 m}) / 2. With I-way interleave (I >= 2) neighbouring cells
belong to different words, bits are i.i.d. with rate lambda(1+f) and

    Pfail = sum_{k>=2} C(72,k) q^k (1-q)^(72-k)  ~  C(72,2) (lambda (1+f) T)^2

With no interleave (I = 1), an MCU puts two flips in the same word; a
72-step recursion over (carried MCU parity, flips so far) gives Pfail
exactly, with first-order term 71 f lambda T. All terms are summed as
positive products, so no cancellation against 1 at tiny rates.

**Design.** Expected uncorrectable words over D days:
E(T) = N Pfail(T) D / T. Without MCU this is linear in T, so scrubbing faster
always helps. With MCU and I = 1 it has a floor 71 f lambda N D that no scrub
rate removes. The longest admissible T is found by bisection and compared
with the closed form T = E_max / (N C(72,2) (lambda(1+f))^2 D).

Environment numbers (2e-10 quiet, 1e-7 events/cell/day for a worst-week SPE,
f = 0.15) are illustrative assumptions, not part data.

## Verification (`make test`)

- H: 72 distinct odd columns, every row weight 26 + 1.
- All 72 single errors corrected and all 2556 double errors flagged on 64
  random words (no data modified on DUE).
- All 59640 triples: 33568 miscorrect, 26072 DUE. Weight-4 codewords
  A4 = 8392 counted two ways (direct 4-subset enumeration and pairs sharing a
  syndrome / 3), and miscorrected triples = 4 A4 exactly.
- Model: I = 1 recursion equals the binomial at f = 0 to 1e-15; the
  C(72,2)(lambda T)^2 and 71 f lambda T asymptotes hold.
- Monte Carlo of a physical array (Poisson upsets, cell-to-word interleave
  map, real decoder) vs. the exact P0, P1, Pfail for four (I, lambda T, f)
  cases, all within 4.5 sigma.

## Run

    make test
    make run

## Sample output

```
Hsiao (72,64) SECDED scrub design: 33554432 words (256 MiB data), f_MCU = 0.15
environment: quiet 2e-10, SPE 1e-07 events/cell/day; requirement E <= 0.01 uncorrectable words

expected uncorrectable words during a 7-day SPE:
  scrub T    no interleave (I=1)   4-way interleave (I=4)   no MCU (f=0 reference)
  1 min             250.1             5.514e-06                4.169e-06
  10 min            250.1             5.514e-05                4.169e-05
  1 h               250.1             0.0003308                0.0002501
  6 h               250.1              0.001985                 0.001501
  1 d               250.2               0.00794                 0.006004
  7 d               250.2               0.05558                  0.04202
  MCU floor at I=1 (T -> 0): 71 f lam N D = 250.1

longest scrub period meeting E <= 0.01:
  case                I=1            I=4            I=4 closed form
  7-day SPE       infeasible     30.23 h        30.23 h      
  SPE scrub load at I=4: 308 words/s (2.8 kB/s of EDAC read bandwidth)
  15-yr quiet     infeasible     402.3 d        402.3 d      
```

Test summary:

```
matrix: 72 distinct odd columns, every row 26 data ones + 1 check
all 4608 single errors corrected, all 163584 double errors flagged DUE (64 random words)
triples: 33568 miscorrect (56.28 %), 26072 DUE; weight-4 codewords A4 = 8392 (direct = pair count / 3, triples = 4*A4)
model: I=1 recursion = binomial at f=0 (worst rel 1.2e-15); C(72,2)(lamT)^2 and 71 f lamT asymptotes hold
Monte Carlo (8192 words x 400 intervals, real decoder) vs exact model:
   I   lamT     f    P0 model  P0 MC     P1 model  P1 MC     Pfail model  Pfail MC    z
   1  2e-03  0.0   0.86601  0.86611   0.12471  0.12458   0.00928      0.00931  +0.52
   1  2e-03  0.3   0.86544  0.86526   0.08838  0.08833   0.04618      0.04641  +1.99
   4  2e-03  0.3   0.82948  0.82972   0.15528  0.15505   0.01524      0.01523  -0.16
   2  8e-03  0.5   0.42366  0.42388   0.36603  0.36570   0.21031      0.21042  +0.48
all checks passed
```

Interleaving, not scrub rate, is the real design driver: with adjacent cells
in the same word, a 7-day SPE leaves about 250 uncorrectable words however
fast you scrub. With 4-way interleave a once-a-day scrub (2.8 kB/s of EDAC
bandwidth) meets the 0.01-word budget, and the quiet-sun limit is over a year.
