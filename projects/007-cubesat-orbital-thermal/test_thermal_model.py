# Astron Built this AP
# Project 007 | Thermal analysis | Transient orbital thermal model of a nadir-pointing 3U CubeSat
# Author: Aditya (aditya-astron)
# Date: 2026-11-02 IST
"""Verification against closed-form results and independent numerics.

Run: python3 -m pytest -q   (or python3 test_thermal_model.py)
"""

import math

import numpy as np

from thermal_model import (ALBEDO, EARTH_IR, R_EARTH, SIGMA, SOLAR_CONST, Orbit,
                           Spacecraft, backward_euler_step, earth_view_factor,
                           eclipse_fraction_analytic, environment_loads,
                           in_cylindrical_shadow, simulate)


def mc_view_factor(lam, h_ratio, n=2_000_000, seed=1):
    """Cosine-weighted rays from the plate; F = fraction that hit Earth.
    A ray hits the sphere iff its angle from nadir is below asin(1/H)."""
    rng = np.random.default_rng(seed)
    u1, u2 = rng.random(n), rng.random(n)
    sin_t = np.sqrt(u1)
    cos_t = np.sqrt(1 - u1)
    ph = 2 * np.pi * u2
    # Local frame: plate normal along local z.
    d = np.stack([sin_t * np.cos(ph), sin_t * np.sin(ph), cos_t], axis=1)
    nadir = np.array([math.sin(lam), 0.0, math.cos(lam)])
    cos_half = math.cos(math.asin(1.0 / h_ratio))
    return float(np.mean(d @ nadir > cos_half))


def test_view_factor_limits():
    H = (R_EARTH + 500e3) / R_EARTH
    assert math.isclose(earth_view_factor(0.0, H), 1 / H**2, rel_tol=1e-12)
    assert earth_view_factor(math.pi, H) == 0.0
    # Continuity at both edges of the partial-visibility band.
    phi = math.asin(1 / H)
    for edge in (math.pi / 2 - phi, math.pi / 2 + phi):
        lo = earth_view_factor(edge - 1e-7, H)
        hi = earth_view_factor(edge + 1e-7, H)
        assert abs(lo - hi) < 1e-5


def test_view_factor_matches_monte_carlo():
    for alt in (400e3, 500e3, 2000e3):
        H = (R_EARTH + alt) / R_EARTH
        for lam_deg in (0, 30, 60, 80, 90, 100, 120):
            lam = math.radians(lam_deg)
            exact = earth_view_factor(lam, H)
            mc = mc_view_factor(lam, H)
            # MC standard error ~ sqrt(F(1-F)/n) <= 3.6e-4
            assert abs(exact - mc) < 2e-3, (alt, lam_deg, exact, mc)


def test_side_plate_reference_value():
    # Vertical plate (lam = 90 deg): F = (phi - sin(phi)cos(phi)) / pi
    # (Gilmore, eq. for a plate perpendicular to the local horizontal).
    H = (R_EARTH + 500e3) / R_EARTH
    phi = math.asin(1 / H)
    ref = (phi - math.sin(phi) * math.cos(phi)) / math.pi
    assert math.isclose(earth_view_factor(math.pi / 2, H), ref, rel_tol=1e-10)


def test_eclipse_fraction_against_closed_form():
    orb = Orbit(alt=500e3, beta=math.radians(20))
    sun = np.array([math.cos(orb.beta), 0, math.sin(orb.beta)])
    th = np.linspace(0, 2 * math.pi, 200_000, endpoint=False)
    shadow = [in_cylindrical_shadow(np.array([math.cos(t), math.sin(t), 0]), sun, orb.radius)
              for t in th]
    assert abs(np.mean(shadow) - eclipse_fraction_analytic(orb.alt, orb.beta)) < 1e-4
    # Hand value for beta = 0 at 500 km: acos(2573.1/6871)/pi = 0.3777
    assert abs(eclipse_fraction_analytic(500e3, 0.0) - 0.3777) < 2e-4
    # Beta beyond asin(R/r) ~ 68 deg: no eclipse.
    assert eclipse_fraction_analytic(500e3, math.radians(70)) == 0.0


def test_nadir_face_loads_at_subsolar_point():
    """At theta = 0, beta = 0 the sun is straight overhead: zenith face (-Z)
    gets full sun, the nadir face (+Z) gets full albedo + IR."""
    sc, orb = Spacecraft(), Orbit(beta=0.0)
    Q, sunlit = environment_loads(sc, orb, 0.0)
    assert sunlit
    F = 1 / ((R_EARTH + orb.alt) / R_EARTH) ** 2
    q_zenith = sc.area[5] * sc.alpha[5] * SOLAR_CONST
    q_nadir = sc.area[4] * (sc.alpha[4] * ALBEDO * SOLAR_CONST * F + sc.eps[4] * EARTH_IR * F)
    assert math.isclose(Q[5], q_zenith, rel_tol=1e-9)
    assert math.isclose(Q[4], q_nadir, rel_tol=1e-9)
    # Eclipse midpoint: no solar, no albedo, IR only.
    Q_ecl, sunlit = environment_loads(sc, orb, math.pi)
    assert not sunlit
    assert Q_ecl[5] == 0.0
    assert math.isclose(Q_ecl[4], sc.area[4] * sc.eps[4] * EARTH_IR * F, rel_tol=1e-9)


def _single_node(T0, dt, n):
    cap, emit = np.array([500.0]), np.array([0.9 * SIGMA * 0.05])
    T = np.array([T0])
    for _ in range(n):
        T = backward_euler_step(T, dt, np.zeros(1), cap, np.zeros((1, 1)), np.zeros((1, 1)), emit)
    return T[0], cap[0], emit[0]


def test_radiative_cooling_closed_form_and_first_order_convergence():
    """m c dT/dt = -k T^4  =>  1/T^3 = 1/T0^3 + 3 k t / (m c)."""
    T0, t_end = 350.0, 4000.0
    errs = []
    for n in (100, 200, 400, 800):
        T_num, C, k = _single_node(T0, t_end / n, n)
        T_ref = (1 / T0**3 + 3 * k * t_end / C) ** (-1 / 3)
        errs.append(abs(T_num - T_ref))
    assert errs[-1] < 0.1
    ratios = [errs[i] / errs[i + 1] for i in range(len(errs) - 1)]
    assert all(1.8 < r < 2.2 for r in ratios), ratios   # backward Euler is O(dt)


def test_isothermal_sphere_equilibrium():
    """Black-body-ish sphere in sunlight, no Earth: T = (alpha S / (4 eps sigma))^(1/4).
    For alpha = eps this is 278.6 K at 1 AU."""
    R = 0.1
    A_tot, A_proj = 4 * math.pi * R**2, math.pi * R**2
    alpha = eps = 0.5
    Q = np.array([alpha * SOLAR_CONST * A_proj])
    emit = np.array([eps * SIGMA * A_tot])
    T = np.array([200.0])
    for _ in range(200):
        T = backward_euler_step(T, 1e4, Q, np.array([1.0]), np.zeros((1, 1)), np.zeros((1, 1)), emit)
    T_ref = (alpha * SOLAR_CONST / (4 * eps * SIGMA)) ** 0.25
    assert abs(T[0] - T_ref) < 1e-6
    assert abs(T_ref - 278.3) < 0.1


def test_periodic_orbit_energy_balance_and_symmetry():
    res = simulate(orb=Orbit(beta=0.0), steps_per_orbit=360)
    # Periodic state: stored energy returns, so absorbed + dissipated = emitted.
    assert abs(res["energy_in"] - res["energy_out"]) / res["energy_in"] < 5e-3
    # beta = 0: the orbit plane is a mirror plane, so +Y and -Y must match.
    assert np.max(np.abs(res["T"][:, 2] - res["T"][:, 3])) < 1e-6
    # Every node within the orbit-average radiative bound.
    assert np.all(res["T"] > 150) and np.all(res["T"] < 400)


def test_orbit_average_against_single_node_estimate():
    """With high conductance the spacecraft is nearly isothermal, so its mean
    temperature must match the orbit-averaged single-node balance
    sum(Q_avg) = sigma sum(eps A) T^4 to within the T^4-averaging error."""
    sc = Spacecraft(g_edge=50.0, g_standoff=50.0)
    orb = Orbit(beta=math.radians(20))
    res = simulate(sc=sc, orb=orb, steps_per_orbit=360)
    th = 2 * math.pi * np.arange(1, 361) / 360
    q_avg = np.mean([environment_loads(sc, orb, t)[0].sum() for t in th])
    T_hand = (q_avg / (SIGMA * np.sum(sc.eps * sc.area))) ** 0.25
    T_mean = res["T"].mean()
    assert abs(T_mean - T_hand) < 1.5, (T_mean, T_hand)


def test_time_step_convergence():
    a = simulate(steps_per_orbit=360)["T"]
    b = simulate(steps_per_orbit=1440)["T"]
    assert abs(a.max() - b.max()) < 1.0
    assert abs(a.min() - b.min()) < 1.0


if __name__ == "__main__":
    import sys
    tests = [v for k, v in dict(globals()).items() if k.startswith("test_")]
    for t in tests:
        t()
        print(f"ok  {t.__name__}")
    print(f"{len(tests)} passed")
    sys.exit(0)
