# Astron Built this AP
# Project 030 | Electronics and circuit analysis | Buck converter: exact periodic steady state, DCM boundary and peak-current-mode subharmonic stability
# Author: Aditya (aditya-astron)
# Date: 2026-11-11 IST
"""Checks the exact switched model against independent references."""
from dataclasses import replace

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from buck import (Buck, ccm_analytic, ccm_report, dcm_ratio, pcm_alpha,
                  pcm_fixed_point, pcm_lyapunov, pcm_map, pcm_run, pss_async, pss_sync,
                  ridley_q, vo_avg_async)

CV = Buck()
D5 = 5.0 / 28.0 * (CV.R + CV.ron + CV.rL) / CV.R


def test_pss_matches_brute_force_ode():
    """Integrate the switched ODE from rest with an adaptive solver for
    many LC time constants; it must land on the matrix-exponential PSS."""
    A, B = CV.matrices()
    x = np.zeros(2)
    for _ in range(600):
        for u, tau in ((CV.vin, D5 * CV.T), (0.0, (1 - D5) * CV.T)):
            sol = solve_ivp(lambda t, y: A @ y + B * u, (0, tau), x,
                            method="DOP853", rtol=1e-11, atol=1e-13)
            x = sol.y[:, -1]
    assert np.allclose(x, pss_sync(CV, D5), rtol=1e-6, atol=1e-8)


def test_ccm_averaged_results():
    sim, ana = ccm_report(CV, D5), ccm_analytic(CV, D5)
    assert sim["vo_avg"] == pytest.approx(5.0, rel=1e-4)
    assert sim["di_pp"] == pytest.approx(ana["di"], rel=1e-3)
    # ESR and capacitive ripple are in quadrature-ish: bounded by max and sum
    assert max(ana["dv_c"], ana["dv_esr"]) < sim["dv_pp"] < ana["dv_c"] + ana["dv_esr"]
    # averaged efficiency ignores ripple RMS; error must be second order
    assert sim["eff"] == pytest.approx(ana["eff"], abs=2e-4)


def test_cap_ripple_formula_without_esr():
    c = replace(CV, rc=0.0)
    sim, ana = ccm_report(c, D5), ccm_analytic(c, D5)
    assert sim["dv_pp"] == pytest.approx(ana["dv_c"], rel=1e-3)


def test_energy_balance():
    sim = ccm_report(CV, D5)
    assert sim["pin"] - sim["pout"] == pytest.approx(sim["loss_ohmic"], rel=1e-9)
    # triangle ripple RMS: Irms^2 = Io^2 + dI^2/12
    io, di = sim["il_avg"], sim["di_pp"]
    assert sim["il_rms2"] == pytest.approx(io ** 2 + di ** 2 / 12, rel=1e-5)


def test_dcm_conversion_ratio():
    c0 = replace(CV, rL=0.0, rc=0.0, ron=0.0)
    for R in (20.0, 50.0, 200.0):
        c = replace(c0, R=R)
        x0, tidle = pss_async(c, 0.25)
        assert tidle > 0
        assert vo_avg_async(c, 0.25, x0) / c.vin == pytest.approx(dcm_ratio(c, 0.25), rel=2e-4)


def test_ccm_dcm_boundary():
    c0 = replace(CV, rL=0.0, rc=0.0, ron=0.0)
    D = 0.25
    Rcrit = 2 * c0.L * c0.fs / (1 - D)
    _, t_below = pss_async(replace(c0, R=0.97 * Rcrit), D)
    _, t_above = pss_async(replace(c0, R=1.03 * Rcrit), D)
    assert t_below == 0.0 and t_above > 0.0


# peak current mode, 5 V -> 3.3 V at 1 MHz
M1, M2, T, IC = (5.0 - 3.3) / 2.2e-6, 3.3 / 2.2e-6, 1e-6, 3.0


@pytest.mark.parametrize("ma", [0.0, 0.3 * M2, 0.5 * M2, 0.9 * M2])
def test_pcm_perturbation_ratio(ma):
    istar = pcm_fixed_point(IC, M1, M2, ma, T)
    assert pcm_map(istar, IC, M1, M2, ma, T) == pytest.approx(istar, abs=1e-12)
    eps = 1e-6
    meas = (pcm_map(istar + eps, IC, M1, M2, ma, T) - istar) / eps
    assert meas == pytest.approx(pcm_alpha(M1, M2, ma), rel=1e-6)


def test_unstable_without_slope_comp():
    """D > 0.5 and no ramp: a 0.1 mA nudge grows by |alpha| per cycle until
    the duty clamps, then wanders inside a bounded band with a positive
    Lyapunov exponent (the stiff-output map is chaotic, not period-2)."""
    istar = pcm_fixed_point(IC, M1, M2, 0.0, T)
    seq = pcm_run(istar + 1e-4, IC, M1, M2, 0.0, T, 2000)
    dev = np.abs(seq[:6] - istar)
    assert np.allclose(dev[1:] / dev[:-1], abs(pcm_alpha(M1, M2, 0.0)), rtol=1e-3)
    tail = seq[500:]
    assert np.ptp(tail) > 0.5
    assert tail.min() >= IC - (M1 + M2) * T - 1e-9 and tail.max() <= IC + M1 * T
    lam = pcm_lyapunov(istar + 1e-4, IC, M1, M2, 0.0, T)
    assert 0.0 < lam <= np.log(abs(pcm_alpha(M1, M2, 0.0)))


def test_half_slope_comp_converges():
    ma = 0.5 * M2
    istar = pcm_fixed_point(IC, M1, M2, ma, T)
    seq = pcm_run(istar + 0.2, IC, M1, M2, ma, T, 60)
    assert abs(seq[-1] - istar) < 1e-12


def test_stability_boundary_matches_ridley():
    # |alpha| = 1 and Ridley's Q -> infinity give the same compensating ramp
    ma = 0.5 * (M2 - M1)
    assert pcm_alpha(M1, M2, ma) == pytest.approx(-1.0, abs=1e-12)
    assert np.isinf(ridley_q(M1, M2, ma))
    # deadbeat ramp ma = m2 gives Ridley's Q = 2/pi
    assert pcm_alpha(M1, M2, M2) == 0.0
    assert ridley_q(M1, M2, M2) == pytest.approx(2 / np.pi, rel=1e-12)
