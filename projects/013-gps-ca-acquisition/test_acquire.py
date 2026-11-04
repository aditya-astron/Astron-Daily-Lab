# Astron Built this AP
# Project 013 | Navigation / Signal processing | GPS L1 C/A code generation and FFT acquisition with CFAR detection
# Author: Aditya (aditya-astron)
# Date: 2026-11-04 IST
"""Verification for the C/A code generator and the acquisition detector.

Run with `python3 test_acquire.py` or `pytest`.
"""
import numpy as np

from acquire import (FS, N_MS, T_CODE, acquire, detection_probability,
                     synth_signal, threshold)
from cacode import CODE_LENGTH, ca_code, ca_code_bits, sample_code

# First 10 chips, octal, leading chip = MSB (IS-GPS-200, Table 3-Ia).
FIRST_10_OCTAL = {
    1: 0o1440, 2: 0o1620, 3: 0o1710, 4: 0o1744, 5: 0o1133, 6: 0o1455,
    7: 0o1131, 8: 0o1454, 9: 0o1626, 10: 0o1504, 11: 0o1642, 12: 0o1750,
    13: 0o1764, 14: 0o1772, 15: 0o1775, 16: 0o1776, 17: 0o1156, 18: 0o1467,
    19: 0o1633, 20: 0o1715, 21: 0o1746, 22: 0o1763, 23: 0o1063, 24: 0o1706,
    25: 0o1743, 26: 0o1761, 27: 0o1770, 28: 0o1774, 29: 0o1127, 30: 0o1453,
    31: 0o1625, 32: 0o1712,
}


def test_first_chips_match_icd():
    for prn, ref in FIRST_10_OCTAL.items():
        bits = ca_code_bits(prn)[:10]
        assert int("".join(map(str, bits)), 2) == ref, prn


def test_gold_correlation_is_three_valued():
    # Periodic auto- and cross-correlation of a 10-stage Gold family take
    # only the values {-65, -1, 63} away from the autocorrelation peak.
    codes = np.array([ca_code(p) for p in range(1, 33)], dtype=float)
    spectra = np.fft.fft(codes, axis=1)
    allowed = {-65, -1, 63}
    for i in range(32):
        corr = np.rint(np.fft.ifft(spectra[i] * np.conj(spectra), axis=1).real).astype(int)
        assert corr[i, 0] == CODE_LENGTH
        assert set(np.unique(corr[i, 1:])) <= allowed, i + 1
        others = np.delete(corr, i, axis=0)
        assert set(np.unique(others)) <= allowed, i + 1


def test_k1_threshold_hand_calc():
    # K = 1: Z ~ Exp(1), so thr = -ln(Pfa_cell) exactly.
    m = 21 * N_MS
    pfa_cell = 1 - (1 - 1e-3) ** (1 / m)
    assert abs(threshold(1e-3, m, 1) / -np.log(pfa_cell) - 1) < 1e-6


def _true_cell_stat(cn0, k, fd, chips, rng):
    """Normalised non-coherent statistic at the exact signal cell, known s2 = 1."""
    sat = {"prn": 7, "cn0_dbhz": cn0, "doppler_hz": fd, "code_phase_chips": chips,
           "phase_rad": rng.uniform(0, 2 * np.pi)}
    r = synth_signal([sat], k, rng=rng)
    t = np.arange(k * N_MS) / FS
    rep = sample_code(7, FS, k * N_MS, -chips) * np.exp(-2j * np.pi * fd * t)
    y = (r * rep).reshape(k, N_MS).sum(axis=1)
    return np.sum(np.abs(y) ** 2) / N_MS


def test_false_alarm_rate_noise_only():
    rng = np.random.default_rng(1)
    k, pfa_cell, trials = 4, 0.05, 4000
    thr = threshold(pfa_cell, 1, k)
    hits = sum(_true_cell_stat(-200.0, k, 0.0, 0.0, rng) > thr for _ in range(trials))
    sigma = np.sqrt(pfa_cell * (1 - pfa_cell) / trials)
    assert abs(hits / trials - pfa_cell) < 4 * sigma, hits / trials


def test_detection_probability_matches_noncentral_chi2():
    rng = np.random.default_rng(2)
    k, trials = 5, 1500
    thr = threshold(1e-3, 21 * N_MS, k)
    for cn0 in (34.0, 36.0):
        pd = detection_probability(cn0, k, thr)
        hits = sum(_true_cell_stat(cn0, k, 1500.0, 333.0, rng) > thr for _ in range(trials))
        sigma = np.sqrt(pd * (1 - pd) / trials)
        assert abs(hits / trials - pd) < 4 * sigma, (cn0, hits / trials, pd)


def test_acquisition_end_to_end():
    rng = np.random.default_rng(42)
    sats = [
        {"prn": 5, "cn0_dbhz": 45.0, "doppler_hz": -4120.0, "code_phase_chips": 250.5},
        {"prn": 19, "cn0_dbhz": 41.0, "doppler_hz": 1875.0, "code_phase_chips": 999.25},
        {"prn": 30, "cn0_dbhz": 38.0, "doppler_hz": 310.0, "code_phase_chips": 3.0},
    ]
    r = synth_signal(sats, 10, rng=rng)
    truth = {s["prn"]: s for s in sats}
    for prn in range(1, 33):
        res = acquire(r, prn, 10)
        if prn not in truth:
            assert not res["detected"], f"false alarm on PRN {prn}"
            continue
        s = truth[prn]
        assert res["detected"], prn
        err = (res["code_phase_chips"] - s["code_phase_chips"] + 511.5) % CODE_LENGTH - 511.5
        assert abs(err) <= 0.125, (prn, err)  # half a sample at 4 samples/chip
        assert abs(res["doppler_hz"] - s["doppler_hz"]) < 15.0, (prn, res["doppler_hz"])
        assert abs(res["cn0_dbhz"] - s["cn0_dbhz"]) < 1.5, (prn, res["cn0_dbhz"])


if __name__ == "__main__":
    tests = [v for n, v in sorted(globals().items()) if n.startswith("test_")]
    for fn in tests:
        fn()
        print(f"PASS {fn.__name__}")
    print(f"{len(tests)} tests passed")
