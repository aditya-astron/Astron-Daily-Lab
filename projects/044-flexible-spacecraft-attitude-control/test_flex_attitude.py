# Astron Built this AP
# Project 044 | PID and control systems | Attitude control of a spacecraft with flexible solar arrays
# Author: Aditya (aditya-astron)
# Date: 2026-11-17 IST
"""Verification. Run: python3 -m pytest -q test_flex_attitude.py"""

import numpy as np
import pytest
from scipy.integrate import quad, solve_ivp

import control as c
import flex_attitude as fa
from flexsat import BETA1, FlexSat, cantilever_mode1

BETA2 = 4.694091132974175


@pytest.fixture(scope="module")
def sat():
    return FlexSat()


def test_cantilever_frequencies(sat):
    # Euler-Bernoulli clamped-free: w_n = beta_n^2 sqrt(EI / (m l^4))
    ref = np.sqrt(sat.EI / (sat.m * sat.l**4))
    assert sat.w_c[0] == pytest.approx(BETA1**2 * ref, rel=1e-7)
    assert sat.w_c[1] == pytest.approx(BETA2**2 * ref, rel=1e-6)


def test_participation_matches_exact_mode_shape(sat):
    norm = quad(lambda x: sat.m * cantilever_mode1(x, sat.l)**2, 0, sat.l)[0]
    L1 = quad(lambda x: sat.m * (sat.r + x) * cantilever_mode1(x, sat.l), 0, sat.l)[0]
    L1 = np.sqrt(2.0) * L1 / np.sqrt(norm)
    assert sat.L[0] == pytest.approx(L1, rel=1e-6)


def test_modal_effective_inertia_is_complete(sat):
    # Parseval in the mass inner product: sum L_k^2 -> panel inertia about the hub axis
    assert np.sum(sat.L**2) == pytest.approx(2 * sat.J_panel, rel=1e-3)


def test_hybrid_model_reproduces_free_free_fe(sat):
    n = len(sat.w_c)
    assert np.allclose(sat.modal_free_free(n)[:3], sat.free_free_fe(3), rtol=1e-9)
    one = sat.w_c[0] / np.sqrt(1 - sat.L[0]**2 / sat.J_t)
    assert sat.modal_free_free(1)[0] == pytest.approx(one, rel=1e-12)
    # three modes are converged on the first frequency to 1e-5
    assert sat.modal_free_free(3)[0] == pytest.approx(sat.free_free_fe(1)[0], rel=1e-5)


def test_zoh_discretisation(sat):
    A, B, _ = sat.state_space(2)
    T = 0.5
    Ad, Bd = c.c2d_zoh(A, B, T)
    x0 = np.linspace(0.1, -0.2, A.shape[0])
    sol = solve_ivp(lambda t, x: A @ x + B[:, 0] * 0.3, (0, T), x0, rtol=1e-11, atol=1e-13)
    assert np.allclose(sol.y[:, -1], Ad @ x0 + Bd[:, 0] * 0.3, atol=1e-9)


def test_rigid_loop_tends_to_continuous_design(sat):
    lp = c.Loop(sat, 1e-3, fa.WN, fa.ZETA, n_modes=0)
    w, z = lp.mode_table()[0]
    assert w == pytest.approx(fa.WN, rel=2e-3)
    assert z == pytest.approx(fa.ZETA, rel=2e-3)


def test_filter_state_space_matches_biquad():
    T = 0.5
    b = np.array([0.3, -0.1, 0.25])
    a = np.array([1.0, -0.9, 0.35])
    f = c.tf2ss_biquad(b, a, "test")
    z = np.exp(1j * np.array([0.0, 0.05, 0.7, 2.0, 3.3, 6.0]) * T)
    ref = np.polyval(b, z) / np.polyval(a, z)      # b0 + b1 z^-1 + b2 z^-2 over a(z^-1)
    assert np.allclose(f.freq(z), ref, rtol=1e-12)


def test_prewarping_puts_features_exactly_where_asked():
    T = 0.5
    n = c.notch(2.5, T, 0.7)
    assert abs(n.freq(np.array([np.exp(2.5j * T)]))[0]) < 1e-12
    lp = c.lowpass2(0.8, T, 0.5)
    H = lp.freq(np.array([np.exp(0.8j * T)]))[0]
    assert abs(H) == pytest.approx(1 / (2 * 0.5), rel=1e-12)
    assert np.degrees(np.angle(H)) == pytest.approx(-90.0, abs=1e-9)


def test_closed_loop_matches_continuous_simulation(sat):
    """Run the plant in continuous time between samples with the held
    command and the flight-code ordering; must equal the augmented matrix."""
    T, wp = 0.5, sat.free_free_fe(1)[0]
    filt = c.notch(0.8 * wp, T, 0.7)
    lp = c.Loop(sat, T, fa.WN, fa.ZETA, filt, n_modes=3)
    X = lp.simulate(fa.STEP, 40)
    A, B, C = lp.A, lp.B, lp.C
    x = np.zeros(A.shape[0])
    xf = np.zeros(filt.A.shape[0])
    u_held = 0.0
    for k in range(40):
        y = C @ x
        e = lp.Kp * (fa.STEP - y[0]) - lp.Kd * y[1]
        u_next = (filt.C @ xf).item() + filt.D * e
        xf = filt.A @ xf + filt.B[:, 0] * e
        sol = solve_ivp(lambda t, s: A @ s + B[:, 0] * u_held, (0, T), x, rtol=1e-11, atol=1e-14)
        x = sol.y[:, -1]
        u_held = u_next
        assert np.allclose(x, X[k + 1, :len(x)], atol=1e-10)
        assert u_held == pytest.approx(X[k + 1, -1], abs=1e-10)


def test_step_has_no_steady_state_error(sat):
    T = 0.5
    lp = c.Loop(sat, T, fa.WN, fa.ZETA, c.lowpass2(0.8, T, 0.5), n_modes=3)
    Phi, G = lp.closed_loop()
    xs = np.linalg.solve(np.eye(len(G)) - Phi, G * fa.STEP)
    assert xs[0] == pytest.approx(fa.STEP, rel=1e-10)
    assert abs(xs[-1]) < 1e-12        # no torque needed at rest


@pytest.mark.parametrize("idx", [1, 2, 3])
def test_nyquist_gain_margin_equals_eigenvalue_crossing(sat, idx):
    name, T, filt = fa.designs(sat)[idx]
    pc, _ = fa.loop(sat, T, filt).margins()
    g_nyq = min(g for _, g in pc if g > 1.0)
    g_eig = fa.gain_band(sat, T, filt)[1]
    assert g_eig == pytest.approx(g_nyq, rel=2e-3), name


def test_design_conclusions(sat):
    d = fa.designs(sat)
    radii = [fa.loop(sat, T, f).spectral_radius() for _, T, f in d]
    assert radii[0] > 1.0 and all(r < 1.0 for r in radii[1:])
    # a notch only phase-stabilises a mode sitting above its centre
    band_on = fa.freq_band(sat, *d[1][1:])
    band_below = fa.freq_band(sat, *d[2][1:])
    assert 0.9 < band_on[0] < 1.0 and band_below[0] < 0.8
    # the plain-PD instability is the first mode alone, not an aliased one
    assert fa.critical_period(sat, 1.0, n_modes=1) < 0.5


def test_conclusions_hold_with_more_modes(sat):
    _, T, filt = fa.designs(sat)[2]
    lo3 = fa.freq_band(sat, T, filt)[0]
    f4 = lambda s: c.Loop(sat, T, fa.WN, fa.ZETA, filt, n_modes=5, freq_scale=s).spectral_radius() - 1
    lo5 = fa.edge(f4, 1.0, -0.01, 0.3)
    assert lo5 == pytest.approx(lo3, abs=0.01)


def test_phase_rule_close_to_exact_at_nominal(sat):
    Tc = fa.critical_period(sat, 1.0)
    Te = fa.phase_rule_period(sat.free_free_fe(1)[0])
    assert abs(Te / Tc - 1) < 0.05
