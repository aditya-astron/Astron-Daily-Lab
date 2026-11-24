# Astron Built this AP
# Project 062 | Signal processing / Embedded firmware | Sinc^N decimation of a second-order sigma-delta stream
# Author: Aditya (aditya-astron)
# Date: 2026-11-24 IST
"""System-level checks of the decimation study (run `make` first)."""
import numpy as np
import pytest

import decimation_study as ds


def test_c_decimator_matches_numpy_boxcar_cascade():
    N, R = 3, 40
    t = np.arange(R * 500) / ds.FS
    bits, _ = ds.modulate(0.4 * np.sin(2 * np.pi * 3e3 * t))
    h = np.ones(1)
    for _ in range(N):
        h = np.convolve(h, np.ones(R))
    ref = np.convolve(bits.astype(np.float64), h)[R - 1:bits.size:R]
    assert np.array_equal(ds.cic(bits, N, R), ref)


def test_output_splits_exactly_into_signal_and_shaped_error():
    # Y = z^-2 X + (1-z^-1)^2 E, so CIC(y) = CIC(x delayed) + CIC(NTF e)
    N, R = 3, 32
    n = R * 400
    x = 0.5 * np.sin(2 * np.pi * 2e3 * np.arange(n) / ds.FS)
    bits, e = ds.modulate(x)
    h = np.ones(1)
    for _ in range(N):
        h = np.convolve(h, np.ones(R))
    sig = np.convolve(np.concatenate([[0, 0], x[:-2]]), h)[R - 1:n:R]
    ne = np.convolve(np.convolve(e, [1.0, -2.0, 1.0])[:n], h)[R - 1:n:R]
    y = ds.cic(bits, N, R)
    # first two outputs see the modulator's zero initial state (e[-1], e[-2])
    assert np.max(np.abs(y[1:] - sig[1:] - ne[1:])) < 1e-6 * R ** N


@pytest.mark.parametrize("R", [8, 32, 100])
def test_noise_gain_closed_forms(R):
    assert ds.noise_gain(1, R) == pytest.approx(4.0)
    assert ds.noise_gain(2, R) == pytest.approx(6.0)
    assert ds.noise_gain(3, R) == pytest.approx(6.0 * R)


def test_snr_slopes_per_octave():
    ratios = [32, 64, 128, 256]
    expect = {1: 6.02, 2: 12.04, 3: 15.05}
    for N, slope in expect.items():
        snr = [ds.run_tone(N, R)[0] for R in ratios]
        fit = np.polyfit(np.log2(ratios), snr, 1)[0]
        assert abs(fit - slope) < 0.5, (N, fit)


@pytest.mark.parametrize("N", [2, 3, 4])
def test_white_error_model_within_1p5_db(N):
    snr, pred, _ = ds.run_tone(N, 64)
    assert abs(snr - pred) < 1.5


def test_quantizer_error_variance_exceeds_uniform_model():
    # e = y - s2 is not uniform on [-1, 1]: the 1-bit loop lets |s2| exceed 1
    _, _, var_e = ds.run_tone(3, 64)
    assert 0.45 < var_e < 0.65


def test_sinc_order_three_meets_requirement_at_chosen_ratio():
    N, R = 3, 81
    assert ds.run_tone(N, R)[0] >= ds.SNR_REQ_DB
    # and N = 2 at the same rate falls well short (aliased shaped noise)
    assert ds.run_tone(2, R)[0] < ds.SNR_REQ_DB - 15


def test_sinc3_bandwidth_and_droop():
    # |H| = |sin(pi f R / fs) / (R sin(pi f / fs))|^3: -11.8 dB, i.e. (2/pi)^3,
    # at the output Nyquist frequency; for large R the -3 dB point tends to
    # the root of sinc(u)^3 = 10^(-3/20), u = pi f R / fs
    R = 81
    fn = ds.FS / (2 * R)
    assert 20 * np.log10(ds.droop(fn, 3, R)) == pytest.approx(
        60 * np.log10(1 / (R * np.sin(np.pi / (2 * R)))), abs=1e-9)
    assert 20 * np.log10(ds.droop(fn, 3, R)) == pytest.approx(60 * np.log10(2 / np.pi), abs=0.01)
    from scipy.optimize import brentq
    u = brentq(lambda u: (np.sin(u) / u) ** 3 - 10 ** (-3 / 20), 0.1, 3.0)
    assert ds.f3db(3, R) * R / ds.FS == pytest.approx(u / np.pi, abs=1e-3)


def test_pruned_18_bit_design():
    N, R = 3, 81
    w, prune = ds.hogenauer_prune(N, R, 18)
    assert w == 22 and prune[0] == 0
    assert all(b2 >= b1 for b1, b2 in zip(prune, prune[1:]))
    base = ds.run_tone(N, R)[0]
    snr = ds.run_tone(N, R, prune=prune)[0]
    pred_loss = 10 * np.log10(1 + ds.pruning_noise_db(N, R, prune)
                              / (0.5 * (ds.AMP * ds.droop(ds.F0, N, R)) ** 2 / 10 ** (base / 10)))
    assert abs((base - snr) - pred_loss) < 0.1
    assert snr >= ds.SNR_REQ_DB
