# 013 - GPS L1 C/A code generation and FFT acquisition with CFAR detection

## Problem

A GPS receiver starting cold does not know which satellites are up, where
each 1023-chip code period starts, or the Doppler shift (up to about
+/-5 kHz for a slow user). Acquisition has to search that space, decide
which satellites are actually present at a controlled false-alarm rate, and
hand the tracking loops a code phase, Doppler and C/N0 good enough to pull in.

## Approach

**Codes.** IS-GPS-200 Gold codes from two 10-stage LFSRs,
G1 = 1 + x^3 + x^10 and G2 = 1 + x^2 + x^3 + x^6 + x^8 + x^9 + x^10, with
the per-PRN G2 phase-selector taps.

**Search.** Complex baseband at 4.092 MHz (4 samples/chip). For each 500 Hz
Doppler bin, wipe the carrier and do the circular correlation over all 4092
code phases with one FFT pair:

    Y(tau) = IFFT( FFT(r e^{-j 2 pi f n / fs}) * conj(FFT(c)) )

**Detector.** Cell power normalised by the noise, Z = |Y|^2 / (N s2), summed
non-coherently over K = 10 ms. Under H0, 2 sum Z ~ chi2(2K); under H1 it is
noncentral with

    lambda = 2 K (C/N0) T

The threshold gives a whole-grid Pfa of 1e-3 via
Pfa_cell = 1 - (1 - Pfa)^(1/M). The noise level is estimated from the median
of the grid (CFAR), so no prior knowledge of s2 is needed.

**Refinement.** A 25 Hz fine Doppler search at the winning code phase, then
the phase slope of consecutive 1 ms prompt correlations,
df = angle(sum P[k] P*[k-1]) / (2 pi T). C/N0 from E[Z] = 1 + (C/N0) T per ms.

**Limits.** No navigation data bit transitions inside the record, code
Doppler neglected (0.03 chip over 10 ms at 5 kHz), integer-sample code
delays in the tests, one front-end noise model.

## Verification

`test_acquire.py` checks:

- first 10 chips of all 32 PRNs against the ICD octal table
- auto/cross-correlation of all 32 codes at every lag is in {-65, -1, 63}
- K = 1 threshold equals the hand result -ln(Pfa_cell)
- Monte Carlo false-alarm rate on noise only (4000 trials, within 4 sigma)
- Monte Carlo Pd at 34 and 36 dB-Hz against the noncentral chi-square
- end to end: three satellites found, no false alarms on the other 29 PRNs,
  code phase within half a sample, Doppler within 15 Hz, C/N0 within 1.5 dB

## Run

    pip install numpy scipy
    python3 test_acquire.py      # or: pytest
    python3 acquire.py

## Sample output

```
fs = 4.092 MHz, 10 ms non-coherent, grid Pfa = 1e-3, Doppler +/-5 kHz @ 500 Hz
PRN    peak    thr   code ph     true     fd Hz    true   C/N0  true
  3   392.3   38.6    101.25   101.25    2210.1  2210.0   46.6  47.0
 11   178.4   38.6    512.00   512.00   -3482.9 -3480.0   42.3  43.0
 17   112.8   38.6    873.75   873.75     641.3   640.0   40.6  40.0
 24    65.2   38.6      7.50     7.50   -1154.2 -1155.0   38.1  37.0

Pd on-grid (no straddle loss), K = 10, grid Pfa = 1e-3:
  32 dB-Hz: Pd = 0.034
  34 dB-Hz: Pd = 0.310
  36 dB-Hz: Pd = 0.887
  38 dB-Hz: Pd = 1.000
  40 dB-Hz: Pd = 1.000
```
