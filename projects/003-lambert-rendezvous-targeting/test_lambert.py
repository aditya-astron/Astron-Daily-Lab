# Astron Built this AP
# Project 003 | Orbital mechanics | Lambert targeting for a LEO rendezvous transfer
# Author: Aditya (aditya-astron)
# Date: 2026-11-01 IST
"""Checks against textbook examples, an independent propagator and Hohmann."""

import numpy as np
import pytest

from lambert import lambert, propagate, stumpff_c, stumpff_s
from rendezvous import Scenario, search


def test_curtis_example_5_2():
    # Curtis, Orbital Mechanics for Engineering Students, Example 5.2 (mu = 398600).
    v1, v2 = lambert([5000, 10000, 2100], [-14600, 2500, 7000], 3600.0, mu=398600.0)
    np.testing.assert_allclose(v1, [-5.9925, 1.9254, 3.2456], atol=1e-3)
    np.testing.assert_allclose(v2, [-3.3125, -4.1966, -0.38529], atol=1e-3)


def test_vallado_example_7_5():
    # Vallado, Fundamentals of Astrodynamics, Example 7-5 (short way, 76 min).
    v1, v2 = lambert([15945.34, 0, 0], [12214.83399, 10249.46731, 0], 76 * 60.0)
    np.testing.assert_allclose(v1, [2.058913, 2.915965, 0], atol=1e-4)
    np.testing.assert_allclose(v2, [-3.451565, 0.910315, 0], atol=1e-4)


@pytest.mark.parametrize("tof,prograde", [
    (600.0, True),      # fast hyperbolic
    (3000.0, True),     # ellipse, short way
    (9000.0, True),     # ellipse, long time
    (5000.0, False),    # long way round (dth > 180)
])
def test_propagated_endpoint_matches(tof, prograde):
    r1 = np.array([7000.0, 500.0, -300.0])
    r2 = np.array([-2000.0, 8500.0, 1200.0])
    v1, v2 = lambert(r1, r2, tof, prograde=prograde)
    r_end, v_end = propagate(r1, v1, tof)
    assert np.linalg.norm(r_end - r2) < 1e-3          # < 1 m
    assert np.linalg.norm(v_end - v2) < 1e-7


def test_energy_and_momentum_consistent():
    r1 = np.array([7000.0, 0, 0])
    r2 = np.array([0, 9000.0, 1000.0])
    v1, v2 = lambert(r1, r2, 2400.0)
    mu = 398600.4418
    e1 = v1 @ v1 / 2 - mu / np.linalg.norm(r1)
    e2 = v2 @ v2 / 2 - mu / np.linalg.norm(r2)
    assert abs(e1 - e2) < 1e-9
    np.testing.assert_allclose(np.cross(r1, v1), np.cross(r2, v2), rtol=1e-11)


def test_stumpff_continuous_through_zero():
    # Both sides of the series/closed-form switch must agree to ~1 ulp.
    for f in (stumpff_c, stumpff_s):
        for zs in (0.1, -0.1):
            assert abs(f(zs * (1 + 1e-12)) - f(zs * (1 - 1e-12))) < 1e-14
    # Spot values: C(pi^2) = 2/pi^2, S(pi^2) = 1/pi^2.
    assert stumpff_c(np.pi**2) == pytest.approx(2 / np.pi**2, rel=1e-14)
    assert stumpff_s(np.pi**2) == pytest.approx(1 / np.pi**2, rel=1e-14)
    assert stumpff_c(0.0) == 0.5 and stumpff_s(0.0) == pytest.approx(1 / 6)


def test_rejects_180_degree_transfer():
    with pytest.raises(ValueError):
        lambert([7000, 0, 0], [-8000, 0, 0], 3000.0)


def test_hohmann_numbers():
    # 300 -> 800 km: hand calc with vis-viva gives 138.16 + 135.68 m/s.
    dv, dv1, dv2, t_h, _ = Scenario().hohmann()
    assert dv1 * 1e3 == pytest.approx(138.16, abs=0.01)
    assert dv2 * 1e3 == pytest.approx(135.68, abs=0.01)
    assert t_h / 60 == pytest.approx(47.82, abs=0.01)


def test_lambert_optimum_recovers_hohmann():
    sc = Scenario()
    dvh, _, _, t_h, tw_h = sc.hohmann()
    tw, tof, dv, *_ = search(sc, n_wait=91, n_tof=61)
    # Hohmann is the coplanar two-impulse optimum for r2/r1 < 11.94,
    # so Lambert can only approach it from above.
    assert dv >= dvh - 1e-9
    assert dv - dvh < 5e-4                  # within 0.5 m/s
    assert abs(tof - t_h) < 60.0
    assert abs(tw - tw_h) < 60.0
