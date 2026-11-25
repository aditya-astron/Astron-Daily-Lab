# Astron Built this AP
# Project 063 | Thermal analysis | Axially grooved aluminium/ammonia heat pipe transport limits
# Author: Aditya (aditya-astron)
# Date: 2026-11-25 IST
"""Checks against reference property data, the exact duct solution, a
finite-difference solve, hand calculation and an ODE march. Run `pytest -q`."""
import numpy as np
import pytest

import heatpipe as hp

C0 = 273.15


@pytest.mark.parametrize("T, key, ref, tol", [
    (0.0, "P", 429.38e3, 2e-3),       # NIST webbook saturation table
    (20.0, "P", 857.48e3, 2e-3),
    (-33.33, "P", 101.325e3, 3e-3),   # normal boiling point
    (0.0, "hfg", 1262.2e3, 3e-3),
    (20.0, "rho_l", 610.2, 2e-3),
    (20.0, "sigma", 0.0212, 0.03),    # correlations differ by ~2 %
])
def test_ammonia_properties(T, key, ref, tol):
    assert getattr(hp.ammonia(T + C0), key) == pytest.approx(ref, rel=tol)


def test_series_matches_finite_difference():
    for w, d in [(0.5e-3, 1.2e-3), (1.0e-3, 1.0e-3), (1.4e-3, 0.4e-3)]:
        assert hp.groove_flow_coeff_fd(w, d, 60) == pytest.approx(
            hp.groove_flow_coeff(w, d), rel=2e-6)


# Shah & London tabulated exact values, and their curve fit (0.06 % worst).
@pytest.mark.parametrize("alpha, table", [
    (1.0, 14.227), (0.5, 15.548), (0.25, 18.233), (0.1, 21.169)])
def test_poiseuille_number_vs_shah_london(alpha, table):
    W, H = alpha, 1.0
    A, Dh = W * H, 2 * W * H / (W + H)
    fre = A * Dh**2 / (2 * hp.rect_flow_coeff(W, H))
    assert fre == pytest.approx(table, abs=1e-3)
    assert fre == pytest.approx(hp.shah_london_fre(alpha), rel=7e-4)


def test_deep_groove_is_a_slot():
    # d >> w: Q/G*mu = w^3 d / 12 * (1 - 0.630 * (w/2)/d) per groove
    w, d = 1e-4, 5e-3
    ref = w**3 * d / 12 * (1 - 0.6302 * (w / 2) / d)
    assert hp.groove_flow_coeff(w, d) == pytest.approx(ref, rel=1e-5)


def test_capillary_limit_hand_calc():
    # Baseline pipe at 20 C with NIST table values typed in by hand.
    rho_l, rho_v, mu_l, mu_v = 610.2, 6.70, 1.33e-4, 9.7e-6
    sigma, hfg = 0.0212, 1186.3e3
    p = hp.Pipe()
    w, d, N = 0.5e-3, 1.2e-3, 31
    assert p.N == N
    C = hp.groove_flow_coeff(w, d)
    dpcap = 2 * sigma / w                                 # 84.8 Pa
    Rl = mu_l / (rho_l * hfg * N * C)                     # Pa/(W m)
    Rv = 128 * mu_v / (np.pi * rho_v * hfg * 8.5e-3**4)
    q_hand = dpcap / ((Rl + Rv) * 0.825)
    assert q_hand == pytest.approx(185.1, rel=2e-3)
    assert p.q_capillary(hp.ammonia(20 + C0)) == pytest.approx(q_hand, rel=0.03)


def test_march_reproduces_effective_length():
    p = hp.Pipe()
    s = hp.ammonia(20 + C0)
    Q = 100.0
    x, Qx, dP = hp.march(p, s, hp.uniform_profile(p, Q), Q)
    R = p.liquid_resistance(s) + p.vapour_resistance(s, Q)
    assert dP[0] == pytest.approx(R * Q * p.Leff, rel=1e-8)
    assert abs(Qx[-1]) < 1e-6
    assert np.argmax(dP) == 0          # worst meniscus at the evaporator end cap


def test_ground_test_line():
    # Laminar 1-g test: Q_max linear in tilt, zero at the static height.
    p = hp.Pipe(w=0.55e-3, d=1.6e-3)
    s = hp.ammonia(20 + C0)
    q = [p.q_capillary(s, tilt_mm=t, one_g=True) for t in (0.0, 1.0, 2.0)]
    assert q[0] - q[1] == pytest.approx(q[1] - q[2], rel=1e-9)
    h_static = (2 * s.sigma / (s.rho_l * hp.G0 * p.w) - p.Dv) * 1e3
    assert p.q_capillary(s, tilt_mm=h_static - 1e-6, one_g=True) < 1e-3
    assert p.q_capillary(s, tilt_mm=h_static + 0.1, one_g=True) == 0.0


@pytest.mark.parametrize("Tc", [-40, -20, 0, 20, 40])
@pytest.mark.parametrize("w", [0.4e-3, 0.8e-3, 1.4e-3])
def test_capillary_limit_is_the_largest_feasible_load(Tc, w):
    p = hp.Pipe(w=w, d=1.6e-3)
    s = hp.ammonia(Tc + C0)
    Q = p.q_capillary(s)
    need = lambda q: q * (p.liquid_resistance(s) + p.vapour_resistance(s, q)) * p.Leff
    drive = p.cap_pressure(s)
    assert need(Q) <= drive * (1 + 1e-9)
    assert need(Q * 1.001) > drive


def test_other_limits_far_above_capillary_in_band():
    p = hp.Pipe(w=0.55e-3, d=1.6e-3)
    for Tc in range(-20, 61, 10):
        L = p.limits(Tc + C0)
        assert min(L, key=L.get) == "capillary_0g"
