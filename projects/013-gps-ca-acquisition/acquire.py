# Astron Built this AP
# Project 013 | Navigation / Signal processing | GPS L1 C/A code generation and FFT acquisition with CFAR detection
# Author: Aditya (aditya-astron)
# Date: 2026-11-04 IST
"""Parallel code-phase acquisition of GPS L1 C/A signals.

Complex baseband model, fs = 4.092 MHz (4 samples/chip), T = 1 ms:

    r[n] = A c[n - tau] exp(j(2 pi fd n / fs + theta)) + w[n],  E|w|^2 = s2

For each Doppler bin the carrier is wiped off and the circular correlation
with the local code over all N code phases is done at once with FFTs:

    Y(tau) = IFFT( FFT(r e^{-j 2 pi f n/fs}) * conj(FFT(c)) )

Normalised cell power Z = |Y|^2 / (N s2). Without signal Z ~ Exp(1); summing
K ms non-coherently, 2 sum Z is chi-square with 2K dof (central under H0,
noncentral with lambda = 2 K (C/N0) T under H1, since the coherent SNR of
one ms is A^2 N / s2 = (C/N0) T with N0 = s2 / fs). The threshold for a
whole-grid false-alarm probability Pfa over M independent-ish cells is set
from the central chi-square with per-cell Pfa_cell = 1 - (1 - Pfa)^(1/M).

s2 is not known to a receiver, so it is estimated from the grid itself
(median of all cells, scaled by the Gamma(K,1) median), which makes this
a CFAR detector.
"""
import numpy as np
from scipy import stats

from cacode import CHIP_RATE, CODE_LENGTH, sample_code

FS = 4.092e6
T_CODE = 1e-3
N_MS = int(round(FS * T_CODE))  # 4092 samples per code period


def synth_signal(sats, n_ms, fs=FS, noise_var=1.0, rng=None):
    """Sum of C/A signals in complex white noise.

    sats: list of dicts with prn, cn0_dbhz, doppler_hz, code_phase_chips,
    phase_rad. Navigation data bits are omitted (no transitions in n_ms).
    Code Doppler is ignored: at 5 kHz it slips 0.03 chip in 10 ms.
    """
    rng = np.random.default_rng() if rng is None else rng
    n = n_ms * int(round(fs * T_CODE))
    t = np.arange(n) / fs
    w = np.sqrt(noise_var / 2) * (rng.standard_normal(n) + 1j * rng.standard_normal(n))
    r = w
    n0 = noise_var / fs
    for s in sats:
        amp = np.sqrt(10 ** (s["cn0_dbhz"] / 10) * n0)
        # Received code at time t started `code_phase` chips earlier, so the
        # local replica sampled at offset -code_phase lines up with it.
        code = sample_code(s["prn"], fs, n, -s["code_phase_chips"])
        r = r + amp * code * np.exp(1j * (2 * np.pi * s["doppler_hz"] * t + s.get("phase_rad", 0.0)))
    return r


def correlation_grid(r, prn, dopplers, n_ms, fs=FS):
    """Non-coherent sum over n_ms of |Y|^2, shape (len(dopplers), N_MS)."""
    n = int(round(fs * T_CODE))
    code_fft = np.conj(np.fft.fft(sample_code(prn, fs, n)))
    t = np.arange(n * n_ms) / fs
    grid = np.zeros((len(dopplers), n))
    for i, f in enumerate(dopplers):
        wiped = (r[: n * n_ms] * np.exp(-2j * np.pi * f * t)).reshape(n_ms, n)
        y = np.fft.ifft(np.fft.fft(wiped, axis=1) * code_fft, axis=1)
        grid[i] = np.sum(np.abs(y) ** 2, axis=0)
    return grid


def threshold(pfa_grid, n_cells, k):
    """Threshold on sum-of-K normalised cell powers for a grid-wide Pfa."""
    pfa_cell = -np.expm1(np.log1p(-pfa_grid) / n_cells)
    return stats.chi2.isf(pfa_cell, 2 * k) / 2


def detection_probability(cn0_dbhz, k, thr, t=T_CODE):
    """Pd for a signal sitting exactly on a grid cell (no straddle loss)."""
    lam = 2 * k * 10 ** (cn0_dbhz / 10) * t
    return stats.ncx2.sf(2 * thr, 2 * k, lam)


def acquire(r, prn, n_ms, pfa=1e-3, fd_max=5000.0, fd_step=500.0, fs=FS):
    """Coarse search plus fine Doppler / C/N0 refinement for one PRN.

    Returns a dict, with detected=False and no estimates if below threshold.
    """
    n = int(round(fs * T_CODE))
    dopplers = np.arange(-fd_max, fd_max + fd_step / 2, fd_step)
    grid = correlation_grid(r, prn, dopplers, n_ms, fs)

    # CFAR noise estimate. Each cell is (N s2) * Gamma(n_ms, 1) under H0.
    noise_scale = np.median(grid) / stats.gamma.median(n_ms)
    z = grid / noise_scale
    thr = threshold(pfa, z.size, n_ms)
    i_f, i_tau = np.unravel_index(np.argmax(z), z.shape)
    result = {"prn": prn, "peak": z[i_f, i_tau], "threshold": thr,
              "detected": bool(z[i_f, i_tau] > thr)}
    if not result["detected"]:
        return result

    # Fine Doppler at the winning code phase: bin width 25 Hz. One-ms
    # coherent loss is sinc^2(df T), so a 250 Hz coarse miss costs 0.9 dB.
    fine = dopplers[i_f] + np.arange(-fd_step, fd_step + 1e-9, 25.0)
    t = np.arange(n * n_ms) / fs
    code = sample_code(prn, fs, n * n_ms, -i_tau * CHIP_RATE / fs)
    zf = np.empty(len(fine))
    for j, f in enumerate(fine):
        prod = (r[: n * n_ms] * code * np.exp(-2j * np.pi * f * t)).reshape(n_ms, n)
        zf[j] = np.sum(np.abs(prod.sum(axis=1)) ** 2) / noise_scale
    j = int(np.argmax(zf))
    fd = fine[j]
    if 0 < j < len(fine) - 1:  # parabolic interpolation of the peak
        a, b, c = zf[j - 1], zf[j], zf[j + 1]
        fd += 0.5 * (a - c) / (a - 2 * b + c) * 25.0

    # Phase-slope refinement: residual carrier rotates the 1 ms prompt
    # correlations by 2 pi df T per block. The angle of sum(P[k] conj P[k-1])
    # is the ML-style frequency estimate (Kay), unambiguous to +/-500 Hz.
    prod = (r[: n * n_ms] * code * np.exp(-2j * np.pi * fd * t)).reshape(n_ms, n)
    prompt = prod.sum(axis=1)
    if n_ms > 1:
        fd += np.angle(np.sum(prompt[1:] * np.conj(prompt[:-1]))) / (2 * np.pi * T_CODE)

    # E[Z] per ms = 1 + (C/N0) T under H1, so (C/N0) T = mean(Z) - 1.
    snr_t = max(zf[j] / n_ms - 1.0, 1e-12)
    result.update({
        "code_phase_chips": (i_tau * CHIP_RATE / fs) % CODE_LENGTH,
        "doppler_hz": fd,
        "cn0_dbhz": 10 * np.log10(snr_t / T_CODE),
    })
    return result


def main():
    rng = np.random.default_rng(13)
    sats = [
        {"prn": 3, "cn0_dbhz": 47.0, "doppler_hz": 2210.0, "code_phase_chips": 101.25, "phase_rad": 0.4},
        {"prn": 11, "cn0_dbhz": 43.0, "doppler_hz": -3480.0, "code_phase_chips": 512.0, "phase_rad": 2.1},
        {"prn": 17, "cn0_dbhz": 40.0, "doppler_hz": 640.0, "code_phase_chips": 873.75, "phase_rad": -1.0},
        {"prn": 24, "cn0_dbhz": 37.0, "doppler_hz": -1155.0, "code_phase_chips": 7.5, "phase_rad": 1.3},
    ]
    n_ms = 10
    r = synth_signal(sats, n_ms, rng=rng)
    truth = {s["prn"]: s for s in sats}

    print(f"fs = {FS/1e6:.3f} MHz, {n_ms} ms non-coherent, grid Pfa = 1e-3, Doppler +/-5 kHz @ 500 Hz")
    print(f"{'PRN':>3} {'peak':>7} {'thr':>6}  {'code ph':>8} {'true':>8}  {'fd Hz':>8} {'true':>7}  "
          f"{'C/N0':>5} {'true':>5}")
    for prn in range(1, 33):
        res = acquire(r, prn, n_ms)
        if not res["detected"] and prn not in truth:
            continue
        s = truth.get(prn)
        line = f"{prn:3d} {res['peak']:7.1f} {res['threshold']:6.1f}"
        if res["detected"]:
            line += (f"  {res['code_phase_chips']:8.2f} {s['code_phase_chips'] if s else float('nan'):8.2f}"
                     f"  {res['doppler_hz']:8.1f} {s['doppler_hz'] if s else float('nan'):7.1f}"
                     f"  {res['cn0_dbhz']:5.1f} {s['cn0_dbhz'] if s else float('nan'):5.1f}")
        else:
            line += "  missed"
        print(line)

    thr = threshold(1e-3, 21 * N_MS, n_ms)
    print("\nPd on-grid (no straddle loss), K = 10, grid Pfa = 1e-3:")
    for cn0 in (32, 34, 36, 38, 40):
        print(f"  {cn0} dB-Hz: Pd = {detection_probability(cn0, n_ms, thr):.3f}")


if __name__ == "__main__":
    main()
