# 010 - Photon-Counting PPM Link for a Deep-Space Optical Downlink

**Domain:** Laser/optical communication and link budgets

## Problem

A spacecraft on a Mars trajectory carries a 4 W, 1550 nm laser behind a 22 cm
telescope. A 5 m ground telescope with a superconducting nanowire detector
counts photons. How many bits per second can the link carry from 0.2 AU out to
2.7 AU, and which pulse-position modulation (PPM) order should it use at each
range?

## Approach

**Link budget (Friis form, optical):**

```
G_tx = (pi D_tx / lambda)^2        G_rx = 4 pi A_rx / lambda^2
P_rx = P_tx G_tx G_rx (lambda / 4 pi R)^2 L_optics L_point L_atm
n_s  = eta P_rx / (h c / lambda)                      [photons/s]
n_b  = eta L_sky A_rx Omega_fov d_lambda L_rx / (h c / lambda)
```

**Poisson M-PPM channel.** Each symbol is M slots of width Ts with one pulsed
slot. At fixed average power the pulse collects `Ks = n_s M Ts` photons; every
slot also sees `Kb = n_b Ts` background photons.

Exact uncoded symbol error rate for max-count detection with random tie break,
summing over the signal-slot count k:

```
Pc = sum_k Pois(k; Ks+Kb) [Fb(k)^M - Fb(k-1)^M] / (M pb(k))
```

(collapsed from the tie-counting binomial sum). With Kb = 0 this reduces to
the erasure result `Pe = (1 - 1/M) e^-Ks`.

Capacity per symbol with uniform inputs, using the slot likelihood ratio
`(1 + Ks/Kb)^k`:

```
C = log2 M - E[ log2 sum_j exp(a (k_j - k_0)) ],   a = ln(1 + Ks/Kb)
```

estimated by Monte Carlo. Throughput is `C / (M Ts)`; Ks is de-rated by 1 dB
as a stand-in for the gap between capacity and a practical SCPPM code. Dead
time, jitter and synchronisation loss are not modelled, so near-range rates
are optimistic (a real detector array saturates first).

## Verification (`test_ppm_link.py`, 23 tests)

- Friis form equals the `P G_tx A / 4 pi R^2` irradiance form; gains, photon
  energy and background rate match hand numbers.
- SER reduces to the erasure closed form as Kb -> 0, equals `1 - 1/M` with no
  signal, and matches a 400k-symbol Monte Carlo within 4 sigma.
- 2-PPM Monte Carlo capacity matches an exact double sum over count pairs.
- M-PPM capacity converges to `(1 - e^-Ks) log2 M` as Kb -> 0, respects
  `0 < C < log2 M`, and satisfies the Fano lower bound from the exact SER.
- Optimal PPM order grows as the link becomes photon-starved.

## Run

```
pip install numpy scipy pytest
python3 ppm_link.py
python3 -m pytest -q
```

## Sample output

```
Link budget at 0.2 AU, lambda = 1550 nm
  Transmit gain [dBi]                  112.98
  Free-space loss [dB]                -347.70
  Receive gain [dBi]                   138.18
  Received power [dBW]                -100.01
  Detected signal photons [1/s]     3.889e+08
  Detected background photons [1/s] 1.930e+07

Best PPM order per range (Ts = 0.5 ns, 1 dB margin to capacity)
  range[AU]  M     Ks      Kb      bits/photon  throughput
       0.2    8    1.236  0.0096        1.30    507.05 Mb/s
       0.5   32    0.791  0.0096        2.28    141.92 Mb/s
       1.0  128    0.791  0.0096        2.85     44.40 Mb/s
       1.5  256    0.703  0.0096        3.04     21.05 Mb/s
       2.0  256    0.395  0.0096        2.90     11.29 Mb/s
       2.7  256    0.217  0.0096        2.51      5.36 Mb/s
```

Photon efficiency peaks near 3 bits/photon around 1.5 AU; beyond that the
order is capped at 256 and background takes a growing share, so efficiency
falls again.
