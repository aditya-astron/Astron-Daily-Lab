# Astron Built this AP
# Project 025 | Signal processing | Polyphase DFT channelizer for a multi-satellite ground receiver
# Author: Aditya (aditya-astron)
# Date: 2026-11-09 IST
"""Verification of the polyphase channelizer. Run: python3 test_channelizer.py (or pytest)."""
import numpy as np
from scipy import signal

import channelizer as ch

FS, M = 8e6, 16
SP = FS / M
H, D_OUT = ch.design_prototype(FS, M)


def resp_db(f):
    _, Hf = signal.freqz(H, worN=[abs(f)], fs=FS)
    return 20 * np.log10(abs(Hf[0]))


def tone_amp(y, f, fs_out, skip=100):
    """Complex amplitude of exp(j 2 pi f t) in y by correlation (steady state)."""
    m = np.arange(len(y))[skip:]
    return np.mean(y[skip:] * np.exp(-2j * np.pi * f * m / fs_out))


def test_polyphase_equals_direct():
    rng = np.random.default_rng(1)
    x = rng.standard_normal(M * 400) + 1j * rng.standard_normal(M * 400)
    y_poly = ch.channelize(x, H, M)
    y_ref = ch.channelize_direct(x, H, M)
    err = np.max(np.abs(y_poly - y_ref)) / np.max(np.abs(y_ref))
    # residual is the rounding of exp(-j 2 pi k n / M) at large n in the reference
    assert err < 1e-10, err


def test_prototype_meets_spec():
    w, Hf = signal.freqz(H, worN=1 << 16, fs=FS)
    mag = 20 * np.log10(np.abs(Hf) + 1e-300)
    assert -mag[w >= 0.6 * SP].max() >= 80.0
    # Kaiser gives equal ripple delta = 10^(-A/20) in pass and stop band:
    # passband peak-to-peak <= 20 log10((1+d)/(1-d)) ~ 0.0017 dB at 80 dB
    d = 10 ** (-80 / 20)
    assert np.ptp(mag[w <= 0.4 * SP]) <= 20 * np.log10((1 + d) / (1 - d)) * 1.05
    # Hand Kaiser estimate N ~ (A - 7.95) / (2.285 dw) + 1, dw = 2 pi (100 kHz / 8 MHz)
    n_hand = (80 - 7.95) / (2.285 * 2 * np.pi * 100e3 / FS) + 1        # 402.4
    n_odd = len(np.trim_zeros(H, "b"))
    assert n_hand <= n_odd < n_hand + 2 * M + 2
    assert (n_odd - 1) // 2 == D_OUT * M and len(H) % M == 0


def test_integer_group_delay():
    x = np.zeros(M * 64, dtype=complex)
    x[0] = 1.0
    y0 = ch.channelize(x, H, M)[0]
    # symmetric FIR: the output peak sits exactly D_OUT samples late
    assert np.argmax(np.abs(y0)) == D_OUT
    assert abs(y0[D_OUT] - H[D_OUT * M]) < 1e-15


def test_tone_routing_and_isolation():
    n = np.arange(M * 3000)
    # (channel, offset from centre): band edges, centre, negative-frequency channel
    cases = [(3, 0.0), (3, +195e3), (3, -195e3), (9, +120e3), (15, -60e3)]
    for k, off in cases:
        f_in = k * SP + off
        x = np.exp(2j * np.pi * f_in * n / FS)
        y = ch.channelize(x, H, M)
        a = tone_amp(y[k], off, SP)
        assert abs(20 * np.log10(abs(a)) - resp_db(off)) < 1e-3
        assert abs(abs(a) - 1) < 2e-4                      # passband gain = 1
        for j in range(M):
            if j != k:
                p = np.mean(np.abs(y[j, 100:]) ** 2)
                assert 10 * np.log10(p) < -80.0, (k, off, j)
    assert abs(ch.channel_freq(15, FS, M) + 500e3) < 1e-6
    assert abs(ch.channel_freq(8, FS, M) + 4e6) < 1e-6


def test_alias_rule_at_band_edge():
    """Content past +250 kHz folds to f - 500 kHz with gain |H(f)|.
    In the transition band (280 kHz) it is a real alias; past 300 kHz it is
    at least 80 dB down, which is why links must stay inside +/-200 kHz."""
    n = np.arange(M * 3000)
    k = 6
    for off, worst_db in [(280e3, None), (320e3, -80.0)]:
        x = np.exp(2j * np.pi * (k * SP + off) * n / FS)
        y = ch.channelize(x, H, M)[k]
        a = tone_amp(y, off - SP, SP)
        got = 20 * np.log10(abs(a))
        assert abs(got - resp_db(off)) < 0.05
        if worst_db is not None:
            assert got < worst_db
    assert resp_db(280e3) > -40                         # transition band leaks


def test_multiply_count():
    direct, poly = ch.mults_per_input_sample(len(H), M)
    assert direct == 16 * (4 + 2 * 432)                  # 13888
    assert poly == (2 * 432 + 4 * 8 * 4) / 16            # 62.0


def test_qpsk_links_through_bank():
    sc = ch.scenario()
    evm = {k: e for k, _, _, e in sc["results"]}
    assert all(e < -55.0 for e in evm.values()), evm
    # weak link between two +30 dB neighbours loses < 10 dB vs. an isolated one
    assert evm[5] - evm[2] < 10.0
    assert sc["floor_dbc"] < -80.0


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"PASS {t.__name__}")
    print(f"{len(tests)} tests passed")
