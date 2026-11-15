# 038 - Carrier-Tracking PLL for an S-band LEO Pass

**Domain:** Signal processing / RF

## Problem

A ground receiver has to hold carrier lock on a 2.25 GHz residual-carrier
downlink from a 500 km circular orbit through an overhead pass, horizon to
horizon above 10 deg elevation. Overhead is the worst case for dynamics:
52 kHz of Doppler, 807 Hz/s of Doppler rate at closest approach and
10 Hz/s^2 of Doppler jerk. Which loop order, what bandwidth, and down to
what C/N0 does it stay locked without a cycle slip?

## Approach

**Geometry.** Non-rotating Earth, station in the orbit plane, so the range and
its first three derivatives are closed form (`th = n t`):

    rho^2 = R^2 + r^2 - 2 R r cos th
    rho rho'' + rho'^2 = R r n^2 cos th
    rho''' = (-R r n^3 sin th - 3 rho' rho'') / rho

Carrier phase `phi = -2 pi fc rho / c`. At TCA `rho'' = R r n^2 / (r - R)`.

**Loop.** Phase-level model of a digital PLL on T-second integrate-and-dump
correlations: the detector is `atan2` of `exp(j(phi - phihat)) + n`, with
`var(n) = 1/((C/N0) T)`. The filter is a chain of forward-Euler accumulators
(order 2: `wn^2, 2 zeta wn`; order 3: `w0^3, 1.1 w0^2, 2.4 w0`, the
Kaplan-Hegarty set with `B_L = 0.7845 w0`). Because every chain ends in
`w0^n`, the steady-state error to a constant n-th phase derivative is exactly
`phi^(n) / w0^n` in discrete time as well.

**Discrete noise bandwidth.** The update is `x' = A x + G v` with
`G = e1 - A e1`; the actual one-sided bandwidth is `B_L = e1' P e1 / 2T`
with `P = A P A' + G G'` (discrete Lyapunov), i.e. `sum(h^2) / 2T`. At
`B_L T = 0.1` the loops are 21-22 % wider than designed.

**Design rule.** `3 sigma + e_dyn <= 45 deg` with `sigma^2 = B_L/(C/N0)`.
For the continuous loop the optimum is
`B_L^(n+1/2) = (2n/3) g sqrt(C/N0)`, `g = phi^(n)_max (B_L/w0)^n`; for the
discrete loop it is found numerically with the actual `B_L` and a stability
check.

**Cycle slips.** The first-order loop gives an exact nonlinear reference:
Tikhonov phase-error density (`E[cos phi] = I1(rho)/I0(rho)`) and the
Viterbi/Lindsey mean slip time `pi^2 rho I0(rho)^2 / (2 B_L)`, which is
also derived independently from the Fokker-Planck first-passage integral.

## Run

    python3 pll.py        # design + full-pass Monte Carlo, ~70 s
    python3 test_pll.py   # verification, ~25 s

Needs numpy and scipy.

## Verification (`test_pll.py`)

- Analog `B_L` of orders 1-3 from `int |H(jw)|^2` matches the closed forms
  (0.25, 0.5303, 0.7845 w0) to 1e-8.
- Lyapunov `B_L` equals the impulse-response `sum(h^2)/2T` to 1e-9 and tends
  to the analog value as `B_L T -> 0`.
- Stability edges: order 1 at `B_L T = 1/2`, order 2 at
  `w0 T = sqrt6 - sqrt2` (Jury) to 1e-9; order 3 impulse response decays
  2 % below and grows 2 % above the numerically found edge.
- Range derivatives vs central differences (1e-6 relative); TCA Doppler
  rate equals the hand formula, 806.64 Hz/s; pass ends at exactly 10 deg.
- Steady-state error `phi''/wn^2` (order 2) and `phi'''/w0^3` (order 3) to
  1e-6; order 3 returns to zero error on a Doppler ramp.
- Noise-free full pass: peak error within 2 % of the quasi-static jerk
  bound.
- Thermal jitter at 45 dB-Hz within 3 % of `B_L,actual / (C/N0)`
  (200 runs x 30 s).
- Mean slip time: closed form vs quadrature to 1e-8 at rho = 2, 5, 10;
  Monte Carlo (2000 slips) within 1-2 % at rho = 1.5 and 2.
- Tikhonov `E[cos phi]`: 0.6988 simulated vs 0.6978.

## Output

    S-band 2.25 GHz, 500 km circular orbit, overhead pass above 10 deg
      pass length              443.0 s
      max Doppler              52.18 kHz
      max Doppler rate         806.6 Hz/s   (at TCA)
      max Doppler jerk         10.18 Hz/s^2

    Design: 3 sigma + e_dyn <= 45 deg at the optimal B_L
      order  T [ms]  C/N0 thr  B_L [Hz]  B_L act  e_dyn [deg]  CN0*T [dB]
          2       1     34.15     88.97   105.79        10.32        4.15
          3       1     21.16      6.48     6.56         6.49       -8.84
          3      10     21.63      6.29     7.07         7.11        1.63

    Full-pass Monte Carlo at the design B_L (40 passes per point)
      order  T [ms]  C/N0  passes w/ slip  sigma sim  sigma lin  max|e| no noise
          2       1   32.1         32/40        18.46      14.55      10.32
          2       1   34.1          1/40        13.69      11.56      10.32
          2       1   36.1          0/40        10.12       9.18      10.32
          3       1   19.2         40/40            -      16.16       6.49
          3       1   21.2          4/40        19.46      12.84       6.49
          3       1   23.2          0/40        13.96      10.20       6.49
          3      10   19.6         10/40        21.27      15.90       7.11
          3      10   21.6          0/40        16.26      12.63       7.11
          3      10   23.6          0/40        12.13      10.03       7.11

## Takeaways

- A second-order loop has to absorb the 807 Hz/s Doppler rate as a static
  error, which forces `B_L ~ 90 Hz` and a 34 dB-Hz threshold. The third-order
  loop only sees jerk and closes the pass at 6.5 Hz, 13 dB lower.
- The 45 deg rule is close but optimistic when the predetection SNR is
  negative: at T = 1 ms (-8.8 dB) the atan detector's noise inflates sigma by
  50 % over the linear figure and 10 % of passes slip at threshold.
  Integrating 10 ms (+1.6 dB) for the same loop gives no slips at threshold,
  provided the carrier is unmodulated over that window (data on a subcarrier).
- Real stations feed Doppler predicts into the NCO, so the loop only tracks
  the orbit-prediction error; the numbers above are the unaided worst case.
