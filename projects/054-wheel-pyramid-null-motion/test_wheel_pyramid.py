# Astron Built this AP
# Project 054 | Attitude control | Four-wheel pyramid momentum envelope and zero-crossing avoidance
# Author: Aditya (aditya-astron)
# Date: 2026-11-21 IST
import numpy as np
import pytest

from wheel_pyramid import (NULL, best_beta, bias_target, fibonacci_sphere,
                           pinv_radius, pyramid, radius_facets, radius_lp,
                           simulate, sphere_radius, zero_crossings)
import run_study as rs

MAGIC = np.arctan(1 / np.sqrt(2))


@pytest.mark.parametrize("deg", [15.0, 30.0, 35.26, 50.0, 70.0])
def test_null_vector(deg):
    assert np.allclose(pyramid(np.deg2rad(deg)) @ NULL, 0, atol=1e-15)


@pytest.mark.parametrize("deg", [15.0, 30.0, 45.0, 70.0])
def test_axis_closed_forms(deg):
    # z: all four wheels add sin(b). x: x = [1,-1,-1,1] hmax gives 4 cos(b)/sqrt2
    # and adding the null vector only raises the largest |x_i|.
    b = np.deg2rad(deg)
    A = pyramid(b)
    assert radius_lp(A, np.array([0, 0, 1.0])) == pytest.approx(4 * np.sin(b), rel=1e-9)
    assert radius_lp(A, np.array([1.0, 0, 0])) == pytest.approx(2 * np.sqrt(2) * np.cos(b), rel=1e-9)


@pytest.mark.parametrize("deg", [20.0, MAGIC * 180 / np.pi, 55.0])
def test_lp_matches_facets(deg):
    A = pyramid(np.deg2rad(deg))
    rng = np.random.default_rng(7)
    for _ in range(200):
        d = rng.normal(size=3)
        d /= np.linalg.norm(d)
        assert radius_lp(A, d) == pytest.approx(radius_facets(A, d), rel=1e-8)


def test_sphere_is_minimum_over_directions():
    for deg in (25.0, 35.26, 60.0):
        A = pyramid(np.deg2rad(deg))
        r = np.array([radius_facets(A, d) for d in fibonacci_sphere(6000)])
        rs_ = sphere_radius(A)
        assert r.min() >= rs_ * (1 - 1e-12)
        assert r.min() <= rs_ * 1.002


def test_best_elevation_is_cube_diagonal():
    # At tan(b) = 1/sqrt2 the axes are (+-1, +-1, 1)/sqrt3. The facet from
    # wheels (1,1,1) and (-1,1,1) has normal (0,-1,1)/sqrt2; the other two
    # wheels project 2/sqrt6 on it, so the inscribed sphere is 4/sqrt6 hmax.
    b, r = best_beta()
    assert b == pytest.approx(MAGIC, abs=1e-6)
    # the optimum is a kink (two facet families cross), so the radius error
    # is first order in the elevation error
    assert r == pytest.approx(4 / np.sqrt(6), rel=1e-8)
    assert sphere_radius(pyramid(MAGIC)) == pytest.approx(4 / np.sqrt(6), rel=1e-14)
    assert np.allclose(pyramid(MAGIC) @ pyramid(MAGIC).T, 4 / 3 * np.eye(3))


def test_min_norm_loss_along_wheel_axis():
    # A^+ = (3/4) A^T at the cube-diagonal angle, A^+ a1 = [3/4, 1/4, -1/4, 1/4]
    # -> reach 4/3. x = [1, 1, -1, 1] reaches 2 (a2 - a3 + a4 = a1). Ratio 2/3,
    # and no direction does worse.
    A = pyramid(MAGIC)
    a1 = A[:, 0]
    assert pinv_radius(A, a1) == pytest.approx(4 / 3, rel=1e-12)
    assert radius_lp(A, a1) == pytest.approx(2.0, rel=1e-9)
    D = fibonacci_sphere(5000)
    ratio = [pinv_radius(A, d) / radius_facets(A, d) for d in D]
    assert min(ratio) >= 2 / 3 - 1e-12


def test_bias_envelope_is_scaled_envelope():
    A = pyramid(np.deg2rad(40.0))
    hmax, hmin = 1.0, 0.15
    lo = np.where(NULL > 0, hmin, -hmax)
    hi = np.where(NULL > 0, hmax, -hmin)
    scale = (hmax - hmin) / (2 * hmax)
    rng = np.random.default_rng(3)
    for _ in range(100):
        d = rng.normal(size=3)
        d /= np.linalg.norm(d)
        assert radius_lp(A, d, lo, hi) == pytest.approx(scale * radius_facets(A, d), rel=1e-8)


def test_bias_target_centres_margins():
    A = pyramid(MAGIC)
    rng = np.random.default_rng(11)
    for _ in range(200):
        h = rng.normal(size=3)
        h *= 0.9 * 0.45 * sphere_radius(A) * rng.uniform() / np.linalg.norm(h)
        p = np.linalg.pinv(A) @ h
        lam, ok = bias_target(p, 0.1, 1.0)
        v = NULL * (p + lam * NULL)
        assert ok
        assert v.min() >= 0.1 - 1e-12 and v.max() <= 1.0 + 1e-12
        assert (v.min() - 0.1) == pytest.approx(1.0 - v.max(), abs=1e-12)


def test_orbit_run():
    A = pyramid(MAGIC)
    common = (A, rs.body_momentum, 3 * rs.PERIOD, 1.0, rs.HMIN, rs.HMAX,
              rs.TAU_NULL, rs.K_NULL)
    t, W0, H, _ = simulate(*common, use_null=False)
    t, W1, H, feas = simulate(*common, use_null=True)
    s = t > 600.0
    assert feas
    assert zero_crossings(W0[s]) > 10
    assert zero_crossings(W1[s]) == 0
    assert np.abs(W1[s]).min() >= rs.HMIN
    assert np.abs(W1).max() <= rs.HMAX
    # null motion never leaks into the body
    assert np.max(np.abs(W1 @ A.T - H)) < 1e-12
    # redistribution torque stays within its allowance
    dW = np.diff(W1 - W0, axis=0)
    assert np.max(np.abs(dW)) <= rs.TAU_NULL + 1e-12
