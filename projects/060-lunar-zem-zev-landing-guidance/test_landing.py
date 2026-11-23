# Astron Built this AP
# Project 060 | Navigation / Guidance | Lunar terminal descent with energy-optimal ZEM/ZEV guidance
# Author: Aditya (aditya-astron)
# Date: 2026-11-23 IST
"""Checks against hand calculations, closed forms and brute force."""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from landing import (G0, G_MOON, Lander, accel_command, energy, feasible_tgo_window,
                     ideal_peak_accel, optimal_tgo, profile_coeffs, simulate)

R0 = np.array([-2000.0, 400.0, 1500.0])
V0 = np.array([45.0, -5.0, -25.0])


def test_vertical_hand_case():
    # 1000 m straight above the pad, at rest, T = 40 s:
    # a0 = -g - 6h/T^2 = 1.62 - 3.75 = -2.13 (z), a1 = 12h/T^3 = 0.1875
    # J_E = a0^2 T + a0 a1 T^2 + a1^2 T^3/3 = 181.476 - 639 + 750 = 292.476
    r, v = np.array([0, 0, 1000.0]), np.zeros(3)
    a0, a1 = profile_coeffs(r, v, 40.0)
    assert np.allclose(a0, [0, 0, -2.13]) and np.allclose(a1, [0, 0, 0.1875])
    assert abs(energy(r, v, 40.0) - 292.476) < 1e-9


def test_profile_reaches_target_and_energy_integral():
    T = 55.0
    a0, a1 = profile_coeffs(R0, V0, T)
    # Exact kinematics of a linear-in-time acceleration plus gravity.
    rT = R0 + V0 * T + 0.5 * (a0 + G_MOON) * T**2 + a1 * T**3 / 6
    vT = V0 + (a0 + G_MOON) * T + 0.5 * a1 * T**2
    assert np.linalg.norm(rT) < 1e-9 and np.linalg.norm(vT) < 1e-12
    JE = quad(lambda t: np.sum((a0 + a1 * t) ** 2), 0, T, epsabs=1e-12)[0]
    assert abs(JE - energy(R0, V0, T)) < 1e-9 * JE


def test_energy_is_minimal_among_feasible_profiles():
    # Adding any control that itself produces zero net r and v change
    # (orthogonal to span{1, t}) must cost more energy.
    T = 55.0
    a0, a1 = profile_coeffs(R0, V0, T)
    base = energy(R0, V0, T)
    # w(t) = t^2 - T t + T^2/6 integrates to 0 and has zero first moment.
    w = lambda t: t * t - T * t + T * T / 6
    assert abs(quad(w, 0, T)[0]) < 1e-9 and abs(quad(lambda t: t * w(t), 0, T)[0]) < 1e-6
    d = np.array([1e-4, -2e-4, 3e-4])
    J = quad(lambda t: np.sum((a0 + a1 * t + d * w(t)) ** 2), 0, T)[0]
    assert J > base


def test_closed_loop_reproduces_open_loop_line():
    # Re-solving at a later point on the arc gives the same straight line.
    T = 55.0
    a0, a1 = profile_coeffs(R0, V0, T)
    for s in [5.0, 20.0, 40.0]:
        r = R0 + V0 * s + 0.5 * (a0 + G_MOON) * s**2 + a1 * s**3 / 6
        v = V0 + (a0 + G_MOON) * s + 0.5 * a1 * s**2
        assert np.allclose(accel_command(r, v, T - s), a0 + a1 * s, atol=1e-10)


def test_optimal_tgo_matches_brute_force():
    for gamma in [0.1, 0.5, 2.0, 10.0]:
        T = optimal_tgo(R0, V0, gamma)
        cost = lambda t: gamma * t + 0.5 * energy(R0, V0, t)
        bf = minimize_scalar(cost, bounds=(5.0, 500.0), method="bounded",
                             options={"xatol": 1e-9}).x
        assert abs(T - bf) < 1e-5, (gamma, T, bf)


def test_unlimited_engine_lands_exactly():
    T = optimal_tgo(R0, V0, 0.5)
    # Every cycle re-solves from the true state and the final segment is
    # flown on the analytic profile, so even a slow 1 Hz ZOH loop lands to
    # round-off: sampling costs delta-v, not accuracy.
    slow = simulate(R0, V0, T, limits=False, dt_guid=1.0, substeps=20)
    assert np.linalg.norm(slow.r_f) < 1e-9 and np.linalg.norm(slow.v_f) < 1e-9
    res = simulate(R0, V0, T, limits=False, dt_guid=0.01, substeps=4)
    assert np.linalg.norm(res.r_f) < 1e-9 and np.linalg.norm(res.v_f) < 1e-9
    # Integrated dv matches the ideal arc's int |u| dt.
    a0, a1 = profile_coeffs(R0, V0, T)
    dv = quad(lambda t: np.linalg.norm(a0 + a1 * t), 0, T)[0]
    assert abs(res.x[-1, 7] - dv) < 0.01


def test_rocket_equation_consistency():
    L = Lander()
    res = simulate(R0, V0, optimal_tgo(R0, V0, 0.5), L)
    m_pred = L.m0 * np.exp(-res.x[-1, 7] / (L.isp * G0))
    assert abs(m_pred - res.x[-1, 6]) < 1e-6


def test_feasible_window_and_nominal_landing():
    L = Lander()
    Tf = feasible_tgo_window(R0, V0, L)
    assert abs(ideal_peak_accel(R0, V0, Tf) - L.t_max / L.m0) < 1e-9
    assert ideal_peak_accel(R0, V0, Tf - 0.5) > L.t_max / L.m0
    res = simulate(R0, V0, Tf, L)
    assert res.sat == 0
    assert np.linalg.norm(res.r_f) < 0.05 and np.linalg.norm(res.v_f) < 0.1
    # Saturating the engine with a too-short T must visibly degrade the landing.
    bad = simulate(R0, V0, 0.7 * Tf, L)
    assert bad.sat > 0 and np.linalg.norm(bad.v_f) > 5.0


if __name__ == "__main__":
    import sys
    tests = [v for k, v in dict(globals()).items() if k.startswith("test_")]
    for t in tests:
        t()
        print("PASS", t.__name__)
    print(f"{len(tests)} tests passed")
    sys.exit(0)
