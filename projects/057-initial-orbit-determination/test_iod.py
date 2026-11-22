# Astron Built this AP
# Project 057 | Orbital mechanics / Navigation | Initial orbit determination: Gibbs, Herrick-Gibbs and Gauss angles-only
# Author: Aditya (aditya-astron)
# Date: 2026-11-22 IST
"""Checks: textbook examples, an independent ODE integrator, convergence
orders, and exact recovery of synthetic orbits. Run with `pytest -q`."""
import numpy as np
import pytest
from scipy.integrate import solve_ivp

import iod
import run_iod

D = np.deg2rad


def test_gibbs_curtis_5_1():
    r1 = np.array([-294.32, 4265.1, 5986.7])
    r2 = np.array([-1365.5, 3637.6, 6346.8])
    r3 = np.array([-2940.3, 2473.7, 6555.8])
    v = iod.gibbs(r1, r2, r3, 398600.0)
    assert np.allclose(v, [-6.2174, -4.0122, 1.5990], atol=1e-4)


def test_herrick_gibbs_vallado_7_4():
    r1 = np.array([3419.85564, 6019.82602, 2784.60022])
    r2 = np.array([2935.91195, 6326.18324, 2660.59584])
    r3 = np.array([2434.95205, 6597.38674, 2521.52311])
    v = iod.herrick_gibbs(r1, r2, r3, 0.0, 76.48, 153.04)
    assert np.allclose(v, [-6.441557, 3.777625, -1.720642], atol=2e-4)


def test_gauss_curtis_5_11_and_5_12():
    th, ra, de, t = run_iod.curtis_gauss_inputs()
    o = run_iod.gauss_curtis(th, ra, de, t)
    assert o["roots"][0] == pytest.approx(9241.8, abs=2.0)
    # Published inputs carry ~6 digits; half-digit rounding alone moves r2
    # by ~1 km (shown in run_iod), so this is the honest tolerance.
    assert np.linalg.norm(o["r2"] - [5662.1, 6538.0, 3269.0]) < 2.0
    assert np.linalg.norm(o["v2"] - [-3.8856, 5.1214, -2.2433]) < 5e-3


def _ode(r0, v0, dt):
    def rhs(_, y):
        r = y[:3]
        return np.concatenate([y[3:], -iod.MU * r / np.linalg.norm(r) ** 3])
    sol = solve_ivp(rhs, (0, dt), np.concatenate([r0, v0]), method="DOP853",
                    rtol=1e-13, atol=1e-12)
    return sol.y[:3, -1], sol.y[3:, -1]


@pytest.mark.parametrize("e,dt", [(0.001, 5000.0), (0.7, 30000.0), (1.4, 20000.0), (0.3, -8000.0)])
def test_universal_propagator_vs_dop853(e, dt):
    a = 9000.0 if e < 1 else -9000.0
    r0, v0 = iod.coe_to_rv(a, e, D(50), D(20), D(70), D(30))
    r, v = iod.propagate(r0, v0, dt)
    rr, vv = _ode(r0, v0, dt)
    assert np.linalg.norm(r - rr) < 1e-5 * np.linalg.norm(rr)
    assert np.linalg.norm(v - vv) < 1e-5 * np.linalg.norm(vv)
    # Lagrange identity f gdot - fdot g = 1 (conservation of h)
    f, g, fd, gd = iod.lagrange_fg(r0, v0, dt)
    assert f * gd - fd * g == pytest.approx(1.0, abs=1e-10)


def test_elements_round_trip():
    el = (7500.0, 0.12, D(63.4), D(210), D(270), D(45))
    back = iod.rv_to_coe(*iod.coe_to_rv(*el))
    assert np.allclose(back, el, atol=1e-10)


def _leo_triplet(sep_deg):
    a = 7078.137
    r0, v0 = iod.coe_to_rv(a, 0.01, D(98.2), D(40), D(60), D(0))
    dt = sep_deg / 360 * 2 * np.pi * np.sqrt(a ** 3 / iod.MU)
    ts = np.array([0.0, dt, 2 * dt])
    rs = [iod.propagate(r0, v0, t)[0] for t in ts]
    return rs, ts, iod.propagate(r0, v0, dt)[1]


def test_herrick_gibbs_truncation_is_fourth_order():
    errs = []
    for sep in (2.0, 4.0, 8.0):
        rs, ts, v2 = _leo_triplet(sep)
        errs.append(np.linalg.norm(iod.herrick_gibbs(*rs, *ts) - v2))
    orders = np.log2(np.array(errs[1:]) / np.array(errs[:-1]))
    assert np.all(np.abs(orders - 4.0) < 0.15)


def test_gibbs_exact_on_any_arc_and_coplanarity_flag():
    for sep in (1.0, 20.0, 60.0):
        rs, _, v2 = _leo_triplet(sep)
        assert np.linalg.norm(iod.gibbs(*rs) - v2) < 1e-9
        assert iod.coplanarity(*rs) < 1e-12
    rs, _, _ = _leo_triplet(10.0)
    n = np.cross(rs[1], rs[2])
    bad = rs[0] + 10.0 * n / np.linalg.norm(n)       # 10 km out of plane
    assert np.degrees(iod.coplanarity(bad, rs[1], rs[2])) > 0.05


def test_gibbs_hg_crossover_inside_classical_band():
    _, cross = run_iod.gibbs_vs_hg(trials=300)
    # Vallado quotes Herrick-Gibbs for separations below about 1-5 deg;
    # with 10 m noise on a 700 km orbit the measured swap is further out.
    assert 5.0 < cross < 20.0


def _leo_pass():
    r0, v0 = iod.coe_to_rv(6378.137 + 800.0, 0.001, D(65.0), D(100), D(30), D(0))
    _, tp, _, _ = run_iod.find_pass(r0, v0)
    return r0, v0, tp


@pytest.mark.parametrize("refine", ["fixed-point", "correct"])
def test_gauss_recovers_exact_leo_orbit(refine):
    r0, v0, tp = _leo_pass()
    ts = [tp - 80, tp, tp + 80]
    L, R = run_iod.observe(r0, v0, ts)
    o = iod.gauss(L, R, ts, refine=refine)
    r2, v2 = iod.propagate(r0, v0, tp)
    assert np.linalg.norm(o["r2"] - r2) < 1e-6        # < 1 mm
    assert np.linalg.norm(o["v2"] - v2) < 1e-8
    # The series start alone is km-level wrong, which is the point of refining.
    assert np.linalg.norm(o["r2_initial"] - r2) > 0.5


def test_long_geo_arc_fixed_point_fails_newton_recovers():
    r0, v0 = iod.coe_to_rv(42164.17, 0.0005, D(3.0), D(80), D(0), D(0))
    tm = 6 * 3600.0
    ts = [tm - 7200, tm, tm + 7200]
    L, R = run_iod.observe(r0, v0, ts)
    with pytest.raises(RuntimeError):
        iod.gauss(L, R, ts, refine="fixed-point")
    o = iod.gauss(L, R, ts, refine="correct")
    r2, v2 = iod.propagate(r0, v0, tm)
    assert np.linalg.norm(o["r2"] - r2) < 1e-6
    assert np.linalg.norm(o["v2"] - v2) < 1e-9


def test_gauss_error_is_linear_in_small_angle_error():
    r0, v0, tp = _leo_pass()
    ts = [tp - 80, tp, tp + 80]
    L, R = run_iod.observe(r0, v0, ts)
    r2 = iod.propagate(r0, v0, tp)[0]
    e1 = np.cross(L[2], [0, 0, 1.0])
    e1 /= np.linalg.norm(e1)
    errs = []
    for k in (1.0, 2.0):
        Lp = [L[0], L[1], L[2] + k * np.radians(1 / 3600) * e1]
        errs.append(np.linalg.norm(iod.gauss(Lp, R, ts, refine="correct")["r2"] - r2))
    assert errs[1] / errs[0] == pytest.approx(2.0, rel=0.02)
