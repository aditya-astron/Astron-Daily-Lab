# Astron Built this AP
# Project 010 | Laser/optical communication and link budgets | Photon-counting PPM link for a deep-space optical downlink
# Author: Aditya (aditya-astron)
# Date: 2026-11-03 IST
"""Checks against hand calculations, closed forms and independent Monte Carlo.

Run: python3 -m pytest -q
"""
import math

import numpy as np
import pytest

from ppm_link import (
    AU,
    OpticalLink,
    best_ppm_order,
    binary_ppm_capacity_exact,
    ppm_capacity,
    ppm_ser,
)


def test_friis_matches_irradiance_form():
    # Friis with G_rx = 4 pi A / lambda^2 must equal P G_tx A / (4 pi R^2).
    link = OpticalLink()
    r = 1.3 * AU
    irradiance_form = link.p_tx * link.tx_gain * link.rx_area / (4 * math.pi * r**2) * link.losses()
    assert link.received_power(r) == pytest.approx(irradiance_form, rel=1e-12)


def test_link_budget_hand_numbers():
    link = OpticalLink()
    # Transmit gain (pi * 0.22 / 1.55e-6)^2 = 112.98 dBi.
    assert 10 * math.log10(link.tx_gain) == pytest.approx(112.98, abs=0.01)
    # 1550 nm photon energy 1.2816e-19 J.
    assert link.photon_energy == pytest.approx(1.2816e-19, rel=1e-4)
    # Inverse-square: doubling range drops the photon rate by exactly 4.
    assert link.signal_photon_rate(AU) / link.signal_photon_rate(2 * AU) == pytest.approx(4.0)
    # Background by hand: 0.1 * 12.566 m^2 * 7.854e-9 sr * 1e-3 um * 0.501 * 0.5 / 1.2816e-19.
    hand = 0.1 * (math.pi * 6.25 * 0.64) * math.pi * 50e-6**2 * 1e-3 * 10 ** -0.3 * 0.5 / 1.2816e-19
    assert link.background_photon_rate() == pytest.approx(hand, rel=1e-3)


@pytest.mark.parametrize("M", [2, 16, 256])
@pytest.mark.parametrize("ks", [0.5, 3.0, 8.0])
def test_ser_noiseless_limit(M, ks):
    # As Kb -> 0 the general formula must approach the erasure result (1 - 1/M) e^-Ks.
    closed = (1 - 1 / M) * math.exp(-ks)
    assert ppm_ser(M, ks, 0.0) == pytest.approx(closed, rel=1e-12)
    assert ppm_ser(M, ks, 1e-7) == pytest.approx(closed, rel=1e-3)


def test_ser_zero_signal_is_random_guess():
    for M in (4, 64):
        assert ppm_ser(M, 0.0, 0.3) == pytest.approx(1 - 1 / M, abs=1e-12)


@pytest.mark.parametrize("M,ks,kb", [(16, 4.0, 0.2), (64, 6.0, 0.05), (8, 2.0, 1.0)])
def test_ser_against_monte_carlo(M, ks, kb):
    rng = np.random.default_rng(7)
    n = 400_000
    counts = rng.poisson(kb, size=(n, M))
    counts[:, 0] += rng.poisson(ks, size=n)
    # Random tie break: add uniform jitter < 1 so it only reorders equal counts.
    decided = np.argmax(counts + rng.random((n, M)), axis=1)
    ser_mc = np.mean(decided != 0)
    sigma = math.sqrt(ser_mc * (1 - ser_mc) / n)
    assert abs(ppm_ser(M, ks, kb) - ser_mc) < 4 * sigma


def test_ser_monotonic():
    sers = [ppm_ser(32, ks, 0.05) for ks in np.linspace(0.5, 12, 12)]
    assert all(a > b for a, b in zip(sers, sers[1:]))


@pytest.mark.parametrize("ks,kb", [(0.5, 0.05), (2.0, 0.2), (4.0, 1.0)])
def test_binary_capacity_mc_matches_exact_sum(ks, kb):
    exact = binary_ppm_capacity_exact(ks, kb)
    mc = ppm_capacity(2, ks, kb, n_mc=400_000, seed=3)
    assert mc == pytest.approx(exact, abs=4e-3)
    assert 0.0 < exact < 1.0


def test_capacity_approaches_erasure_channel():
    for M in (16, 128):
        ks = 1.5
        erasure = (1 - math.exp(-ks)) * math.log2(M)
        assert ppm_capacity(M, ks, 0.0) == pytest.approx(erasure)
        assert ppm_capacity(M, ks, 1e-5, n_mc=100_000) == pytest.approx(erasure, rel=5e-3)


def test_capacity_bounds_and_background_penalty():
    clean = ppm_capacity(64, 3.0, 0.001)
    noisy = ppm_capacity(64, 3.0, 0.2)
    assert 0 < noisy < clean < math.log2(64)


def test_capacity_respects_fano_against_ser():
    # Fano: H(X|Y) <= h(Pe) + Pe log2(M-1), so C >= log2 M - h(Pe) - Pe log2(M-1).
    M, ks, kb = 32, 3.0, 0.05
    pe = ppm_ser(M, ks, kb)
    h = -pe * math.log2(pe) - (1 - pe) * math.log2(1 - pe)
    lower = math.log2(M) - h - pe * math.log2(M - 1)
    assert ppm_capacity(M, ks, kb) >= lower - 0.01


def test_photon_starved_link_prefers_higher_order():
    link = OpticalLink()
    nb = link.background_photon_rate()
    m_near, *_ = best_ppm_order(link.signal_photon_rate(0.2 * AU), nb, 0.5e-9, n_mc=20_000)
    m_far, *_ = best_ppm_order(link.signal_photon_rate(2.7 * AU), nb, 0.5e-9, n_mc=20_000)
    assert m_far > m_near
