# 059 - Optimal ASM frame synchronization: Massey rule vs correlation

## Problem

A CCSDS telemetry receiver finds frame boundaries by locating the 32-bit
attached sync marker 0x1ACFFC1D in the symbol stream. Project 002 did this
with a hard-decision threshold (d_H <= 2), which is fine at uncoded SNRs.
With LDPC or turbo coding the channel symbols arrive at Es/N0 around -2 dB,
the hard BER is above 10 %, and the question becomes which statistic to
search with. This project builds the exact error statistics of the marker
(including windows that partly overlap it), designs the hard threshold, and
compares one-shot acquisition over a full RS(255,223) I=5 CADU
(L = 10232 candidate offsets) for:

| rule | window metric |
|---|---|
| hard | number of hard decisions that agree with the ASM |
| corr | sum s_i y_i |
| massey | sum s_i y_i - (1/lam) sum ln cosh(lam y_i), lam = 2 Es/N0 |
| massey-hiSNR | sum s_i y_i - sum \|y_i\| = -2 sum over sign mismatches of \|y_i\| |
| massey(+-3 dB est) | Massey with the SNR estimate off by a factor of 2 |

## Approach

**Channel.** BPSK, unit-energy symbols, y = x + n with n ~ N(0, N0/2);
randomized frame data modelled as independent equiprobable symbols;
p = Q(sqrt(2 Es/N0)).

**Massey's rule.** Per sample, p(y | ASM symbol s) / p(y | random data) =
exp(lam y s) / cosh(lam y). The window log-likelihood ratio is therefore
lam * [sum s y - (1/lam) sum ln cosh(lam y)]. Correlation alone rewards any
window sitting on large-amplitude samples, whatever their sign; the ln cosh
term removes that bias. As lam -> infinity it becomes sum |y|, giving an
SNR-free rule.

**Overlap windows.** A window shifted k bits from a marker sees 32-k marker
bits, m_k of which disagree with the reference (m_k from
popcount((ASM << k) XOR ASM) on the top 32-k bits), plus k data bits. Its
Hamming distance is exactly

    D_k = Bin(32-k-m_k, p) + Bin(m_k, 1-p) + Bin(k, 1/2)

and its correlation is (32 - 2m_k - 2B) + N(0, 32 sigma^2), B ~ Bin(k, 1/2).
Inside a CADU search there are 31 such windows on each side (the tail of
this marker and the head of the next one).

**Acquisition probability (analytic).** Treat the L window scores as
independent, using the exact pmfs above for the 62 overlap windows. For
correlation: P = int f_T(t) prod_j F_j(t) dt. For hard decisions with random
tie-breaking, given D_T = d and K windows tied at d the win probability is
E[1/(1+K)] = int_0^1 prod_j (P(D_j > d) + P(D_j = d) u) du. Massey has no
convenient closed form and is evaluated by Monte Carlo only.

## Checks (`sync_test.go`)

| Check | Result |
|---|---|
| ASM bit order vs the binary string of 0x1ACFFC1D | exact |
| m_k by bit trick vs character comparison, k = 1..31 | exact; m_1 = 11 bit transitions |
| Random-window false alarm at E = 2 | 529 / 2^32 = 1.232e-7 (project 002 figure) |
| Overlap distance pmf: sum and mean vs (32-k-m)p + m(1-p) + k/2 | 1e-12 |
| Overlap distance pmf vs 200k simulated shifted windows (k = 1, 7, 24) | every bin within 5 sigma |
| Massey metric times lam vs LLR from the Gaussian densities | 1e-9 |
| High-SNR limit: penalty within 32 ln2 / lam of sum \|y\| | holds for lam = 1, 10, 1000 |
| CADU acquisition, hard and corr: Monte Carlo vs analytic at -4 and 0 dB | within 4 sigma |
| Massey fails less than both hard and corr at -4 and 0 dB | yes |

Hand check: at Es/N0 = 0 dB, p = Q(sqrt 2) = 0.0786.

## Results

- The marker's worst aperiodic sidelobe is +9 at a one-bit shift (main peak
  32). A noiseless one-bit slip already needs 11 errors to pass, so for
  E <= 6 the overlap windows add under 0.5 % to the per-frame false-alarm
  count; the random-data windows dominate.
- Hard threshold: E = 2 gives 1.3e-3 false locks per frame search but misses
  81 % of markers at -2 dB. No threshold works at coded-link SNRs; one has to
  search for the best window instead.
- Best-window search over one CADU (5000 trials per point, analytic hard and
  corr in the `th` columns agree with Monte Carlo to within sampling error):

| Es/N0 for P(fail) | hard | corr | massey | massey-hiSNR |
|---|---|---|---|---|
| 0.10 | -0.63 dB | +0.08 dB | -2.55 dB | -2.31 dB |
| 0.01 | +1.40 dB | +2.84 dB | -0.74 dB | -0.62 dB |

- Plain soft correlation is worse than hard decisions above about -2.5 dB:
  among 10^4 random windows some land on unusually large samples, and only
  the correlator is fooled by that. Massey's rule is 2.6 dB better than
  correlation at 10 % failure and 3.6 dB better at 1 %.
- It is robust to the SNR estimate: a factor of 2 either way costs under
  0.2 dB, and the estimate-free -sum|y| form costs about 0.1-0.25 dB.

## Run

    go test ./...                # unit checks + Monte Carlo acquisition check (~2 s)
    go run .                     # full study, ~1 min on 2 cores
    go run . -trials 1000        # quicker

Go 1.22+ (uses math/rand/v2), standard library only. Full output in
`sample_output.txt`.

## Sample output

    Es/N0 for P(fail) = 0.01 (interpolated):
      hard               +1.40 dB   (-1.43 dB vs corr)
      corr               +2.84 dB   (+0.00 dB vs corr)
      massey             -0.74 dB   (-3.58 dB vs corr)
      massey-hiSNR       -0.62 dB   (-3.45 dB vs corr)
      massey(-3dB est)   -0.57 dB   (-3.40 dB vs corr)
      massey(+3dB est)   -0.68 dB   (-3.51 dB vs corr)

Reference: J. L. Massey, "Optimum Frame Synchronization", IEEE Trans.
Communications, COM-20(2), 1972.
