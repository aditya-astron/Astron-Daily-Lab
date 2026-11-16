# Astron Built this AP
# Project 040 | Orbital mechanics | Earth-Moon L2 halo family to the 9:2 NRHO
# Author: Aditya (aditya-astron)
# Date: 2026-11-16 IST
"""Verification for the halo/NRHO pipeline. Run: python3 test_halo_nrho.py"""

import numpy as np

import cr3bp as c
import families as f
import halo_nrho as h
import richardson as r

DAY = 86400.0
_cache = {}


def nrho():
    if "nrho" not in _cache:
        _, _, h0 = h.halo_from_richardson(10000.0)
        _cache["fam"], _cache["nrho"] = h.walk_to_nrho(h0)
    return _cache["nrho"]


def bifurcation():
    if "bif" not in _cache:
        _cache["bif"] = h.halo_bifurcation(h.lyapunov_family())
    return _cache["bif"]


def test_libration_points():
    # Reference values for mu = 0.01215 (Szebehely; Koon et al.)
    g1, g2, x1, x2 = c.libration_points()
    assert abs(x1 - 0.836915) < 2e-6
    assert abs(x2 - 1.155682) < 2e-6
    assert abs(c.potential_gradient(np.array([x1, 0, 0]))).max() < 1e-13
    assert abs(c.potential_gradient(np.array([x2, 0, 0]))).max() < 1e-13
    assert abs(c.jacobi(np.array([x1, 0, 0, 0, 0, 0])) - 3.18834) < 1e-5
    assert abs(c.jacobi(np.array([x2, 0, 0, 0, 0, 0])) - 3.17216) < 1e-5


def test_richardson_sun_earth_table():
    # Richardson (1980) Sun-Earth L1 coefficients, mu = 3.040357e-6
    rc = r.coefficients(1, mu=3.040357143e-6)
    ref = dict(c2=4.0610735, c3=3.0200105, c4=3.0305427, lam=2.0864535,
               k=3.2292681, delta=0.2922144, l1=-15.965599, l2=1.7409005,
               s1=-0.8246605, s2=0.1210986)
    for name, val in ref.items():
        got = getattr(rc, name)
        assert abs(got - val) < 2e-6 * max(1, abs(val)), (name, got, val)
    # ISEE-3 class halo: Az = 120,000 km gives Ax ~ 206,000 km, ~178 day period
    L = 1.49598e8
    Ax = r.ax_from_az(rc, 120000 / (rc.gamma * L)) * rc.gamma * L
    assert abs(Ax - 206000) < 2000
    _, T = r.halo_state(rc, 120000, L, 0.0)
    assert abs(T / (2 * np.pi) * 365.25 - 177.7) < 0.5


def test_richardson_velocity_is_time_derivative():
    rc = r.coefficients(2)
    s, T = r.halo_state(rc, 20000, c.L_STAR, 0.7)
    dtau = 1e-6
    sp, _ = r.halo_state(rc, 20000, c.L_STAR, 0.7 + dtau)
    sm, _ = r.halo_state(rc, 20000, c.L_STAR, 0.7 - dtau)
    rate = 2 * np.pi / T                     # d tau1 / dt
    v_fd = (sp[:3] - sm[:3]) / (2 * dtau) * rate
    assert np.allclose(v_fd, s[3:], atol=1e-9)


def test_stm_matches_finite_differences():
    s0 = np.array([1.08, 0.0, -0.15, 0.01, -0.12, 0.02])
    _, Phi, _ = c.propagate(s0, 1.0, stm=True)
    h_ = 1e-6
    Phi_fd = np.zeros((6, 6))
    for j in range(6):
        e = np.zeros(6)
        e[j] = h_
        p, _ = c.propagate(s0 + e, 1.0)
        m, _ = c.propagate(s0 - e, 1.0)
        Phi_fd[:, j] = (p - m) / (2 * h_)
    assert np.abs(Phi - Phi_fd).max() < 1e-6 * np.abs(Phi).max()


def test_stm_is_symplectic():
    s0 = np.array([1.08, 0.0, -0.15, 0.01, -0.12, 0.02])
    _, Phi, _ = c.propagate(s0, 2.0, stm=True)
    Om = c.symplectic_form()
    assert np.abs(Phi.T @ Om @ Phi - Om).max() < 1e-9 * np.abs(Phi).max() ** 2
    assert abs(np.linalg.det(Phi) - 1) < 1e-8


def test_monodromy_from_half_matches_full_propagation():
    orb = nrho()
    _, Phi_T, _ = c.propagate(orb.s0, orb.period, stm=True)
    M = orb.monodromy()
    assert np.abs(M - Phi_T).max() < 1e-6 * np.abs(Phi_T).max()


def test_lyapunov_bifurcation():
    bif = bifurcation()
    assert abs(f.lyapunov_out_of_plane_trace(bif) - 2.0) < 1e-8
    # Commonly reported Earth-Moon L2 halo bifurcation: C ~ 3.152, T ~ 14.8 d
    assert abs(bif.C - 3.1521) < 2e-4
    assert abs(bif.period * c.T_STAR / DAY - 14.83) < 0.05
    # Third-order theory puts the same branch point at Ax_min
    rc = r.coefficients(2)
    ax_min = np.sqrt(-rc.delta / rc.l1) * rc.gamma * c.L_STAR
    ax_bif = (bif.s_half[0] - bif.s0[0]) / 2 * c.L_STAR
    assert abs(ax_bif / ax_min - 1) < 0.03


def test_halo_family_emanates_from_bifurcation():
    # Halo orbits approach the bifurcating Lyapunov orbit as Az -> 0, with
    # Jacobi constant and period changing quadratically in Az.
    bif = bifurcation()
    dC, dT = [], []
    for Az in (1000.0, 2000.0):
        _, _, orb = h.halo_from_richardson(Az)
        dC.append(bif.C - orb.C)
        dT.append(bif.period - orb.period)
    assert 0 < dC[0] < 1e-4
    assert abs(dC[1] / dC[0] - 4) < 0.1
    assert abs(dT[1] / dT[0] - 4) < 0.1


def test_richardson_guess_close_to_corrected_halo():
    guess, T_guess, orb = h.halo_from_richardson(10000.0)
    assert abs(guess[0] - orb.s0[0]) * c.L_STAR < 1000
    assert abs(T_guess / orb.period - 1) < 0.005
    # Corrected orbit is periodic: back to the start after one period.
    end, _ = c.propagate(orb.s0, orb.period)
    assert np.linalg.norm(end - orb.s0) < 1e-8


def test_nrho_92_geometry_and_stability():
    orb = nrho()
    T_days = orb.period * c.T_STAR / DAY
    assert abs(T_days - 29.530589 / 4.5) < 1e-6
    rp, ra = h.perilune_apolune(orb)
    # Published 9:2 NRHO: rp ~ 3,200-3,500 km, ra ~ 71,000 km,
    # stability index ~ 1.3 (Davis, Zimovan-Spreen et al.)
    assert 3000 < rp < 3600
    assert 70000 < ra < 72500
    assert orb.s_half[2] > 0 and orb.s0[2] < 0     # perilune north, apolune south
    ev, nus = f.stability(orb)
    assert 1.2 < abs(nus[0]) < 1.45
    assert abs(nus[1]) < 1                         # second mode oscillatory
    # Reciprocal pairs: product of all eigenvalues is det(M) = 1
    assert abs(np.prod(ev).real - 1) < 1e-8
    assert np.sort(np.abs(ev - 1))[1] < 1e-3       # trivial pair at +1


def test_nrho_periodic_and_jacobi_conserved():
    orb = nrho()
    end, sol = c.propagate(orb.s0, 3 * orb.period, dense=True)
    C0 = orb.C
    t = np.linspace(0, 3 * orb.period, 3001)
    Cs = np.array([c.jacobi(s) for s in sol.sol(t).T])
    assert np.abs(Cs - C0).max() < 1e-10
    # Mildly unstable (|lambda| ~ 2.2/rev), so three periods still close.
    assert np.linalg.norm(end[:3] - orb.s0[:3]) * c.L_STAR < 1.0     # km


def test_continuation_satisfies_constraints():
    nrho()
    fam = _cache["fam"]
    for q in fam:
        assert abs(q.s_half[3]) < 1e-10 and abs(q.s_half[5]) < 1e-10
        assert abs(q.s_half[1]) < 1e-12
    steps = [np.linalg.norm(fam[i + 1].s0[[0, 2, 4]] - fam[i].s0[[0, 2, 4]])
             for i in range(len(fam) - 2)]
    assert max(steps) < 0.0102 and min(steps) > 0.0098


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"PASS  {t.__name__}")
    print(f"\n{len(tests)} tests passed")
