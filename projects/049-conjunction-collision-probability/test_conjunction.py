# Astron Built this AP
# Project 049 | Orbital mechanics | Conjunction collision probability and avoidance burn sizing
# Author: Aditya (aditya-astron)
# Date: 2026-11-20 IST
"""Checks against closed forms, scipy's noncentral chi-square, Monte Carlo and
two-body propagation."""

import math

import numpy as np
import pytest
from scipy.stats import ncx2

import conjunction as cj
import run_conjunction as rc


def rot(a):
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, -s], [s, c]])


def aniso(sx, sy, ang):
    q = rot(ang)
    return q @ np.diag([sx * sx, sy * sy]) @ q.T


METHODS = [cj.pc_foster, cj.pc_1d, cj.pc_chan]


# -- isotropic closed forms --------------------------------------------------

def test_hand_value_centered_isotropic():
    # sigma = 100 m, R = 20 m: Pc = 1 - exp(-0.02) = 0.0198013266932...
    c = np.eye(2) * 100.0 ** 2
    for f in METHODS:
        assert f(np.zeros(2), c, 20.0) == pytest.approx(0.019801326693244747, rel=1e-12)


@pytest.mark.parametrize("sig, d, r", [
    (100.0, 0.0, 5.0), (100.0, 150.0, 10.0), (50.0, 300.0, 20.0),
    (200.0, 50.0, 300.0), (30.0, 200.0, 10.0),
])
def test_isotropic_offset_is_noncentral_chi2(sig, d, r):
    # |x|^2/sigma^2 is chi'^2 with 2 dof and noncentrality d^2/sigma^2.
    ref = ncx2.cdf((r / sig) ** 2, 2, (d / sig) ** 2)
    m = np.array([d * 0.6, -d * 0.8])
    c = np.eye(2) * sig * sig
    for f in (cj.pc_1d, cj.pc_chan):
        assert f(m, c, r) == pytest.approx(ref, rel=1e-9)
    # Foster at default resolution is meant for R up to a few sigma.
    if r <= 3.0 * sig:
        assert cj.pc_foster(m, c, r) == pytest.approx(ref, rel=1e-9)


def test_small_object_hand_value():
    sig, d, r = 120.0, 90.0, 2.0
    exp = (r * r / (2.0 * sig * sig)) * math.exp(-d * d / (2.0 * sig * sig))
    assert cj.pc_small(np.array([d, 0.0]), np.eye(2) * sig ** 2, r) == pytest.approx(exp, rel=1e-14)


# -- anisotropic agreement ---------------------------------------------------

@pytest.mark.parametrize("seed", range(12))
def test_foster_matches_1d_random(seed):
    rng = np.random.default_rng(seed)
    sy = rng.uniform(20.0, 200.0)
    sx = sy * rng.uniform(1.0, 30.0)
    c = aniso(sx, sy, rng.uniform(0, math.pi))
    m = rng.normal(0.0, 1.5, 2) * np.array([sx, sy])
    r = rng.uniform(1.0, 2.0 * sy)
    pf, p1 = cj.pc_foster(m, c, r), cj.pc_1d(m, c, r)
    assert pf == pytest.approx(p1, rel=1e-9)


def test_foster_converges_exponentially():
    # R = 2.5 sigma_min is a hard case for the default grid.
    c, m, r = aniso(400.0, 40.0, 0.7), np.array([60.0, -30.0]), 100.0
    ref = cj.pc_1d(m, c, r)
    errs = [abs(cj.pc_foster(m, c, r, n_r=n, n_theta=2 * n) / ref - 1) for n in (2, 4, 8, 16)]
    assert errs[-1] < 1e-14
    assert errs[1] < 0.05 * errs[0] and errs[2] < 1e-4 * errs[1]


def test_rotation_invariance():
    c, m, r = aniso(300.0, 50.0, 0.0), np.array([80.0, 40.0]), 12.0
    base = cj.pc_1d(m, c, r)
    for a in (0.3, 1.1, 2.5):
        q = rot(a)
        for f in METHODS:
            assert f(q @ m, q @ c @ q.T, r) == pytest.approx(f(m, c, r), rel=1e-11)
        assert cj.pc_1d(q @ m, q @ c @ q.T, r) == pytest.approx(base, rel=1e-11)


def test_chan_error_grows_with_aspect():
    m, r = np.array([50.0, 20.0]), 15.0
    errs = []
    for ar in (1.0, 3.0, 10.0, 30.0):
        c = aniso(60.0 * ar, 60.0, 0.4)
        errs.append(abs(cj.pc_chan(m, c, r) / cj.pc_1d(m, c, r) - 1.0))
    assert errs[0] < 1e-12
    assert errs[1] < errs[2] < errs[3]


def test_small_object_limit():
    c, m = aniso(250.0, 70.0, 0.9), np.array([100.0, 40.0])
    ratios = [cj.pc_small(m, c, r) / cj.pc_1d(m, c, r) for r in (20.0, 5.0, 1.0, 0.1)]
    assert abs(ratios[-1] - 1.0) < 1e-5
    assert all(abs(a - 1) > abs(b - 1) for a, b in zip(ratios, ratios[1:]))


def test_large_body_saturates():
    c = aniso(50.0, 20.0, 0.2)
    assert cj.pc_1d(np.array([30.0, 10.0]), c, 2000.0) == pytest.approx(1.0, abs=1e-12)


# -- covariance scaling ------------------------------------------------------

def test_max_pc_matches_closed_form_for_small_object():
    c, m, r = aniso(300.0, 60.0, 0.5), np.array([250.0, 90.0]), 0.5
    k_num, p_num = cj.max_pc_numeric(m, c, r)
    k_an, p_an = cj.max_pc_small_object(m, c, r)
    assert k_num == pytest.approx(k_an, rel=1e-4)
    assert p_num == pytest.approx(p_an, rel=1e-4)


# -- encounter geometry ------------------------------------------------------

def test_encounter_frame_and_projection():
    r1, v1, r2, v2, c1, c2 = rc.scenario()
    dr, dv = r2 - r1, v2 - v1
    e = cj.encounter_frame(dr, dv)
    assert np.allclose(e @ e.T, np.eye(3), atol=1e-14)
    assert np.linalg.det(e) == pytest.approx(1.0, abs=1e-14)
    m, c = cj.project_encounter(dr, dv, c1 + c2)
    _, d = cj.rectilinear_miss(dr, dv)
    assert np.linalg.norm(m) == pytest.approx(d, rel=1e-12)
    assert abs(m[1]) < 1e-9 * d                      # miss lies on the x axis
    assert np.all(np.linalg.eigvalsh(c) > 0.0)
    # RTN covariance round trip preserves the trace.
    assert np.trace(c1) == pytest.approx(np.sum(rc.SIG_PRIMARY ** 2), rel=1e-12)


def test_scenario_monte_carlo():
    r1, v1, r2, v2, c1, c2 = rc.scenario()
    dr, dv = r2 - r1, v2 - v1
    cov = c1 + c2
    m, c = cj.project_encounter(dr, dv, cov)
    pc = cj.pc_1d(m, c, rc.HBR)
    n = 1_500_000
    p_mc, _ = cj.pc_monte_carlo(dr, dv, cov, np.eye(3) * 2 * rc.SIG_VEL ** 2,
                                rc.HBR, n, np.random.default_rng(7))
    assert abs(p_mc - pc) < 4.0 * math.sqrt(pc * (1 - pc) / n)


# -- avoidance manoeuvre -----------------------------------------------------

def test_cw_whole_orbit_is_pure_along_track_drift():
    n = 1.05e-3
    t = 3 * 2 * math.pi / n
    d = cj.cw_impulse_displacement(n, t, (0.0, 0.02, 0.0))
    assert d[0] == pytest.approx(0.0, abs=1e-9)
    assert d[1] == pytest.approx(-3.0 * t * 0.02, rel=1e-12)


@pytest.mark.parametrize("lead_orb, dv_rtn", [
    (0.5, (0.0, -0.05, 0.0)), (1.0, (0.0, 0.03, 0.0)),
    (2.0, (0.02, 0.0, 0.0)), (1.25, (0.0, 0.0, 0.05)),
])
def test_cw_matches_two_body(lead_orb, dv_rtn):
    r1, v1, *_ = rc.scenario()
    a = np.linalg.norm(r1)
    n = math.sqrt(cj.MU_EARTH / a ** 3)
    lead = lead_orb * 2 * math.pi / n
    d_cw = cj.cw_impulse_displacement(n, lead, dv_rtn)
    d_tb = cj.twobody_impulse_displacement(r1, v1, lead, dv_rtn)
    # Linearisation error scales like |d|^2 / a.
    assert np.linalg.norm(d_cw - d_tb) < 2e-3 * np.linalg.norm(d_tb) + 1e-3
    assert np.linalg.norm(d_cw - d_tb) < 10.0 * np.linalg.norm(d_tb) ** 2 / a + 1e-3
