# Astron Built this AP
# Project 048 | Electronics and circuit analysis | ADC anti-alias filter: Chebyshev LC ladder vs. Sallen-Key cascade
# Author: Aditya (aditya-astron)
# Date: 2026-11-20 IST
"""Checks against published tables, scipy.signal and independent solves."""

import math

import numpy as np
import pytest
import scipy.signal as ss

import filters as fl
import run_filter as rf

FP = 8e6
WP = fl.TWO_PI * FP
F = np.linspace(1e3, 60e6, 1501)


def ladder(n, ap, g=None):
    g = g or fl.cheby_g(n, ap)
    elems, rl = fl.ladder_values(g, 50.0, FP)
    return [k for k, _ in elems], np.array([v for _, v in elems]), rl


# -- prototype -------------------------------------------------------------

@pytest.mark.parametrize("n, ap, ref", [
    # Matthaei, Young & Jones, Table 4.05-2(a)
    (5, 0.5, [1.7058, 1.2296, 2.5408, 1.2296, 1.7058, 1.0000]),
    (4, 0.5, [1.6703, 1.1926, 2.3661, 0.8419, 1.9841]),
    (3, 0.1, [1.0316, 1.1474, 1.0316, 1.0000]),
    (3, 0.01, [0.6292, 0.9703, 0.6292, 1.0000]),
])
def test_cheby_g_matches_table(n, ap, ref):
    assert np.allclose(fl.cheby_g(n, ap)[1:], ref, atol=6e-5)


def test_butter_g_matches_table():
    assert np.allclose(fl.butter_g(5)[1:-1], [0.6180, 1.6180, 2.0, 1.6180, 0.6180], atol=1e-4)


def test_order_is_minimum():
    n = fl.cheby_order(0.5, 60.0, 4.0)
    assert n == 5
    att = lambda m: -10 * math.log10(fl.cheby_mag2(4.0, m, 0.5))
    assert att(5) >= 60.0 > att(4)
    nb = fl.butter_order(0.5, 60.0, 4.0)
    eps2 = 10 ** 0.05 - 1
    att_b = lambda m: 10 * math.log10(1 + eps2 * 4.0 ** (2 * m))
    assert nb == 6 and att_b(6) >= 60.0 > att_b(5)


def test_poles_match_scipy():
    z, p, k = ss.cheb1ap(5, 0.5)
    assert np.allclose(np.sort_complex(fl.cheby_poles(5, 0.5)), np.sort_complex(p), atol=1e-12)


# -- ladder ----------------------------------------------------------------

@pytest.mark.parametrize("n", [3, 5, 7])
def test_odd_ladder_equals_chebyshev_transfer_function(n):
    kinds, vals, rl = ladder(n, 0.5)
    s21, _ = fl.ladder_sparams(F, kinds, vals, 50.0, rl)
    b, a = ss.cheby1(n, 0.5, WP, analog=True)
    _, h = ss.freqs(b, a, fl.TWO_PI * F)
    assert np.max(np.abs(s21 - h)) < 1e-11            # magnitude and phase
    assert np.allclose(np.abs(s21) ** 2, fl.cheby_mag2(F / FP, n, 0.5), atol=1e-12)


def test_even_ladder_needs_unequal_load_and_matches_scipy():
    kinds, vals, rl = ladder(4, 0.5)
    assert rl == pytest.approx(50.0 / 1.9841, rel=1e-4)   # g5 is a conductance here
    s21, _ = fl.ladder_sparams(F, kinds, vals, 50.0, rl)
    b, a = ss.cheby1(4, 0.5, WP, analog=True)
    _, h = ss.freqs(b, a, fl.TWO_PI * F)
    assert np.max(np.abs(s21 - h)) < 1e-11
    assert 20 * math.log10(abs(s21[0])) == pytest.approx(-0.5, abs=1e-6)


def test_butterworth_ladder_is_3db_at_edge():
    kinds, vals, rl = ladder(5, None, fl.butter_g(5))
    s21, _ = fl.ladder_sparams([FP], kinds, vals, 50.0, rl)
    assert abs(s21[0]) ** 2 == pytest.approx(0.5, abs=1e-12)


def test_lossless_power_conservation_and_reflection_zeros():
    kinds, vals, rl = ladder(5, 0.5)
    s21, s11 = fl.ladder_sparams(F, kinds, vals, 50.0, rl)
    assert np.max(np.abs(np.abs(s21) ** 2 + np.abs(s11) ** 2 - 1)) < 1e-12
    frz = np.cos((2 * np.arange(1, 4) - 1) * math.pi / 10) * FP
    _, s11z = fl.ladder_sparams(frz, kinds, vals, 50.0, rl)
    assert np.max(np.abs(s11z)) < 1e-9


def test_group_delay_from_poles_matches_ladder_phase():
    kinds, vals, rl = ladder(5, 0.5)
    f = np.linspace(0.05e6, 7.95e6, 50)
    df = 1.0
    ph = lambda ff: np.angle(fl.ladder_sparams(ff, kinds, vals, 50.0, rl)[0])
    tau_num = -np.angle(np.exp(1j * (ph(f + df) - ph(f - df)))) / (fl.TWO_PI * 2 * df)
    tau = fl.pole_group_delay(fl.TWO_PI * f, fl.cheby_poles(5, 0.5) * WP)
    assert np.allclose(tau_num, tau, rtol=1e-6)


# -- dissipation -----------------------------------------------------------

def test_uniform_q_loss_is_delay_times_omega_over_q():
    kinds, vals, rl = ladder(5, 0.5)
    q = 1e4
    f = np.array([1e6, 4e6, 7e6, 8e6])
    lossy = fl.ladder_sparams(f, kinds, vals, 50.0, rl, q, q)[0]
    clean = fl.ladder_sparams(f, kinds, vals, 50.0, rl)[0]
    d_db = 20 * np.log10(np.abs(clean / lossy))
    tau = fl.pole_group_delay(fl.TWO_PI * f, fl.cheby_poles(5, 0.5) * WP)
    pred = 20 / math.log(10) * fl.TWO_PI * f * tau / q
    assert np.allclose(d_db, pred, rtol=2e-3)


def test_phase_sensitivity_rule_for_unequal_q():
    kinds, vals, rl = ladder(5, 0.5)
    ql, qc = 5e3, 5e4
    f = np.array([2e6, 6e6, 8e6])
    lossy = fl.ladder_sparams(f, kinds, vals, 50.0, rl, ql, qc)[0]
    clean = fl.ladder_sparams(f, kinds, vals, 50.0, rl)[0]
    sens = rf.ladder_phase_sens(f, kinds, vals, rl)
    qs = np.array([ql if k == "L" else qc for k in kinds])[:, None]
    pred = -20 / math.log(10) * np.sum(sens / qs, axis=0)
    assert np.allclose(20 * np.log10(np.abs(clean / lossy)), pred, rtol=2e-3)
    # sum of phase sensitivities over all reactive elements = -w tau (frequency scaling)
    tau = fl.pole_group_delay(fl.TWO_PI * f, fl.cheby_poles(5, 0.5) * WP)
    assert np.allclose(sens.sum(axis=0), -fl.TWO_PI * f * tau, rtol=1e-6)


# -- Orchard ---------------------------------------------------------------

def test_zero_first_order_sensitivity_at_reflection_zeros():
    kinds, vals, rl = ladder(5, 0.5)
    frz = np.cos((2 * np.arange(1, 3) - 1) * math.pi / 10) * FP
    s = rf.mag_sens_db(lambda x: fl.ladder_sparams(frz, kinds, x, 50.0, rl)[0], vals)
    assert np.max(np.abs(s)) < 1e-8
    # any perturbation can only lower |S21| below 1 (passive, matched source)
    rng = np.random.default_rng(1)
    x = vals * (1 + 0.1 * rng.uniform(-1, 1, (500, 5)))
    assert np.max(np.abs(fl.ladder_sparams(F, kinds, x, 50.0, rl)[0])) <= 1 + 1e-12


# -- Sallen-Key ------------------------------------------------------------

def sk_mna(f, r1, r2, c1, c2, gbw):
    """Independent nodal solve: unknowns Vx, Vp, Vo; op-amp Vo = A (Vp - Vo)."""
    out = []
    for ff in f:
        s = 1j * fl.TWO_PI * ff
        a = fl.TWO_PI * gbw / s
        m = np.array([
            [1 / r1 + 1 / r2 + s * c1, -1 / r2, -s * c1],
            [-1 / r2, 1 / r2 + s * c2, 0],
            [0, -a, 1 + a],
        ])
        out.append(np.linalg.solve(m, np.array([1 / r1, 0, 0]))[2])
    return np.array(out)


def test_sallen_key_cascade_equals_ladder_response():
    _, _, elems, rl, poles, secs = rf.nominal_designs()
    kinds = [k for k, _ in elems]
    vals = np.array([v for _, v in elems])
    s21, _ = fl.ladder_sparams(F, kinds, vals, 50.0, rl)
    assert np.max(np.abs(fl.sallen_key_response(F, secs) - s21)) < 1e-11
    for s in secs:
        if s[0] == "SK":
            q = 0.5 * math.sqrt(s[3] / s[4])
            assert q in (pytest.approx(1.1783, abs=1e-3), pytest.approx(4.5450, abs=1e-3))


def test_finite_gbw_stage_against_mna():
    r1, r2, c1, c2 = 200.0, 220.0, 888e-12, 10.75e-12
    f = np.linspace(1e5, 50e6, 40)
    for gbw in (3e8, 4e9):
        h = fl.sallen_key_response(f, [("SK", r1, r2, c1, c2)], gbw)
        assert np.allclose(h, sk_mna(f, r1, r2, c1, c2, gbw), rtol=1e-10)
    assert np.allclose(fl.sallen_key_response(f, [("SK", r1, r2, c1, c2)], 1e15),
                       fl.sallen_key_response(f, [("SK", r1, r2, c1, c2)]), atol=1e-6)


def test_min_gbw_is_the_threshold():
    secs = rf.nominal_designs()[5]
    g = rf.min_gbw(secs)
    ideal = rf.db(fl.sallen_key_response(rf.F_PB, secs))
    err = lambda gb: np.max(np.abs(rf.db(fl.sallen_key_response(rf.F_PB, secs, gb)) - ideal))
    assert err(g) <= 0.1 < err(0.98 * g)


def test_monte_carlo_sanity():
    rows = dict(rf.monte_carlo(400, seed=3))
    ripple, rej = rows["LC ladder, L +-2 %, C +-2 %"]
    assert np.all(ripple >= 0) and np.all(rej > 70)
    ripple5, _ = rows["LC ladder, L/C +-5 %"]
    assert np.median(ripple5) > np.median(ripple) > 0.5 - 1e-3
