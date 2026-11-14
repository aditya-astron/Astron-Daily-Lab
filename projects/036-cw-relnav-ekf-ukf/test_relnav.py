# Astron Built this AP
# Project 036 | Kalman filtering and sensor fusion | Rendezvous relative navigation in the Hill frame: EKF vs UKF vs the PCRB
# Author: Aditya (aditya-astron)
# Date: 2026-11-14 IST
"""Verification for relnav.py. Run: python3 test_relnav.py"""

import numpy as np
from scipy.integrate import quad_vec, solve_ivp
from scipy.linalg import expm, solve_discrete_are
from scipy.stats import chi2

from relnav import (Scenario, UT, cw_A, cw_phi, van_loan_Q, h_lidar, H_lidar,
                    h_angles, H_angles, predict, ekf_update, ukf_update,
                    rts_smooth, info_bound, monte_carlo, range_sigma,
                    mean_motion)

N_MM = mean_motion(420e3)
H_POS = np.hstack([np.eye(3), np.zeros((3, 3))])


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {detail}")
    assert ok, name


def test_stm():
    A = cw_A(N_MM)
    worst = 0.0
    for t in (1.0, 10.0, 600.0, 2900.0, 5578.0, 20000.0):
        P = cw_phi(N_MM, t)
        worst = max(worst, np.abs(P - expm(A * t)).max() / np.abs(P).max())
    semi = np.abs(cw_phi(N_MM, 700.0) @ cw_phi(N_MM, 1300.0) - cw_phi(N_MM, 2000.0)).max()
    check("CW STM equals expm(A t) and is a semigroup", worst < 1e-12 and semi < 1e-9,
          f"(rel err {worst:.1e}, semigroup {semi:.1e})")
    # Against direct integration of the CW ODE from a drift state.
    x0 = np.array([-500.0, -8000.0, 100.0, 0.02, 0.85, -0.01])
    sol = solve_ivp(lambda t, x: A @ x, (0, 6000.0), x0, rtol=1e-12, atol=1e-12)
    d = np.abs(sol.y[:, -1] - cw_phi(N_MM, 6000.0) @ x0)
    check("CW STM matches numerical integration", d[:3].max() < 1e-6, f"({d[:3].max():.1e} m)")


def test_van_loan():
    q, dt = 1e-10, 10.0
    G = np.vstack([np.zeros((3, 3)), np.eye(3)])
    ref, _ = quad_vec(lambda s: cw_phi(N_MM, s) @ G @ G.T @ cw_phi(N_MM, s).T * q,
                      0.0, dt, epsrel=1e-12)
    err = np.abs(van_loan_Q(N_MM, dt, q) - ref).max() / np.abs(ref).max()
    # Double-integrator limit: Q_vv = q dt, Q_rv = q dt^2/2, Q_rr = q dt^3/3,
    # with CW corrections of order (n dt)^2 = 1.3e-4.
    Q = van_loan_Q(N_MM, dt, q)
    di = np.array([Q[0, 0] / (q * dt**3 / 3), Q[0, 3] / (q * dt**2 / 2), Q[3, 3] / (q * dt)])
    check("Van Loan Q equals quadrature and the double-integrator limit",
          err < 1e-10 and np.all(np.abs(di - 1) < 2 * (N_MM * dt) ** 2), f"(rel err {err:.1e}, DI ratios {np.round(di, 6)})")


def test_jacobian():
    rng = np.random.default_rng(3)
    worst = 0.0
    for _ in range(200):
        x = np.r_[rng.normal(0, 3000, 3) + [0, -5000, 0], rng.normal(0, 1, 3)]
        Hn = np.zeros((3, 6))
        for j in range(3):
            d = np.zeros(6)
            d[j] = 1e-3
            Hn[:, j] = (h_lidar(x + d) - h_lidar(x - d)) / 2e-3
        H = H_lidar(x)
        worst = max(worst, np.abs(H - Hn).max() / np.abs(H).max())
    check("lidar Jacobian matches central differences", worst < 1e-7, f"(rel {worst:.1e})")


def test_angles_only_observability():
    """Angles are degree-0 homogeneous; CW is linear => scale is unobservable
    on a coast. Linearised: H(x_k) Phi_k x0 = H(x_k) x_k = 0 (Euler)."""
    sc = Scenario(lidar=False)
    x0 = sc.x0
    same = max(np.abs(h_angles(cw_phi(sc.n, t) @ (1.7 * x0)) - h_angles(cw_phi(sc.n, t) @ x0)).max()
               for t in np.linspace(0, sc.T, 50))
    Ri = np.linalg.inv(sc.R)

    def gramian(burn):
        W = np.zeros((6, 6))
        x, Phi = x0.copy(), np.eye(6)
        for k in range(sc.N + 1):
            H = H_angles(x)
            W += (H @ Phi).T @ Ri @ (H @ Phi)
            x = sc.Phi @ x + (np.r_[0, 0, 0, 0.1, 0, 0] if burn and k == sc.N // 2 else 0)
            Phi = sc.Phi @ Phi
        return W

    # Scale states to comparable units (y = S^-1 x) before comparing eigenvalues.
    S = np.diag([1, 1, 1, sc.n, sc.n, sc.n])     # velocities as v/n metres
    for burn in (False, True):
        Ws = S @ gramian(burn) @ S
        w = np.linalg.eigvalsh(Ws)
        y = np.linalg.solve(S, x0)               # x0 in scaled coordinates
        ray = np.linalg.norm(Ws @ y) / (w[-1] * np.linalg.norm(y))
        if not burn:
            check("coast: scaled trajectory gives identical angles", same < 1e-12, f"({same:.1e} rad)")
            check("coast: initial state spans the Gramian null space (range scale)",
                  ray < 1e-13 and w[1] / w[-1] > 1e3 * max(ray, 1e-16),
                  f"(|W x0|/(lambda_max |x0|) {ray:.1e}, next eigenvalue ratio {w[1] / w[-1]:.1e})")
        else:
            check("radial burn: Gramian becomes full rank", w[0] / w[-1] > 1e-9,
                  f"(lambda_min/max {w[0] / w[-1]:.1e})")


def test_ukf_equals_kf_linear():
    """UT is exact for linear maps, so on a linear model the UKF is the KF."""
    sc = Scenario()
    rng = np.random.default_rng(11)
    R = np.eye(3) * 4.0
    h = lambda x: H_POS @ x
    ut = UT(6)
    xk, Pk = sc.x0.copy(), sc.P0.copy()
    xu, Pu = sc.x0.copy(), sc.P0.copy()
    worst = 0.0
    for k in range(300):
        z = h(sc.x0) + rng.normal(0, 50, 3)
        xk, Pk = predict(xk, Pk, sc.Phi, sc.Q)
        xu, Pu = predict(xu, Pu, sc.Phi, sc.Q)
        xk, Pk, _, _ = ekf_update(xk, Pk, z, h, lambda x: H_POS, R, [])
        xu, Pu, _, _ = ukf_update(xu, Pu, z, h, R, [], ut)
        worst = max(worst, np.abs(xk - xu).max() / np.abs(xk).max(),
                    np.abs(Pk - Pu).max() / np.abs(Pk).max())
    check("UKF reproduces the linear KF exactly", worst < 1e-9, f"(rel {worst:.1e})")


def test_dare():
    sc = Scenario()
    R = np.eye(3) * 1.0
    P = sc.P0.copy()
    for _ in range(20000):
        P = sc.Phi @ P @ sc.Phi.T + sc.Q
        P = P - P @ H_POS.T @ np.linalg.solve(H_POS @ P @ H_POS.T + R, H_POS @ P)
    Pp = sc.Phi @ P @ sc.Phi.T + sc.Q
    X = solve_discrete_are(sc.Phi.T, H_POS.T, sc.Q, R)
    err = np.abs(Pp - X).max() / np.abs(X).max()
    check("Riccati recursion converges to the DARE solution", err < 1e-8, f"(rel {err:.1e})")


def test_rts_equals_batch():
    """No process noise: the smoothed initial state equals the batch MAP."""
    sc = Scenario()
    rng = np.random.default_rng(5)
    R = np.eye(3) * 25.0
    N = 200
    Q0 = np.zeros((6, 6))
    xt = sc.x0 + rng.multivariate_normal(np.zeros(6), sc.P0)
    zs, x = [], xt.copy()
    for k in range(N):
        zs.append(H_POS @ x + rng.normal(0, 5, 3))
        x = sc.Phi @ x
    xf, Pf, xp, Pp = [], [], [], []
    x, P = sc.x0.copy(), sc.P0.copy()
    for k in range(N):
        if k:
            x, P = predict(x, P, sc.Phi, Q0)
        xp.append(x); Pp.append(P)
        x, P, _, _ = ekf_update(x, P, zs[k], lambda v: H_POS @ v, lambda v: H_POS, R, [])
        xf.append(x); Pf.append(P)
    xs, Ps = rts_smooth(np.array(xf), np.array(Pf), np.array(xp), np.array(Pp), sc.Phi)
    J = np.linalg.inv(sc.P0)
    b = J @ sc.x0
    Phi = np.eye(6)
    for k in range(N):
        M = H_POS @ Phi
        J += M.T @ M / 25.0
        b += M.T @ zs[k] / 25.0
        Phi = sc.Phi @ Phi
    xb = np.linalg.solve(J, b)
    ex = np.abs(xs[0] - xb).max()
    eP = np.abs(Ps[0] - np.linalg.inv(J)).max() / np.abs(Ps[0]).max()
    check("RTS smoothed x0 equals the batch MAP solution", ex < 1e-6 and eP < 1e-6,
          f"(state {ex:.1e}, cov rel {eP:.1e})")


def test_monte_carlo_consistency():
    runs = 50
    sc = Scenario(lidar=True)
    truths, res = monte_carlo(sc, runs, seed=21)
    lo, hi = chi2.ppf([0.005, 0.995], runs * 6) / runs
    lo3, hi3 = chi2.ppf([0.005, 0.995], runs * 3) / runs
    bound = np.mean([np.trace(info_bound(sc, X)[:, :3, :3], axis1=1, axis2=2) for X in truths], axis=0)
    tail = slice(sc.N // 3, None)
    for kind in ("ekf", "ukf"):
        anees = res[kind]["nees"].mean(0)
        anis = res[kind]["nis"].mean(0)
        f1 = np.mean((anees > lo) & (anees < hi))
        f2 = np.mean((anis > lo3) & (anis < hi3))
        mse = (res[kind]["err"][:, :, :3] ** 2).sum(-1).mean(0)
        eff = np.mean(mse[tail] / bound[tail])
        check(f"{kind.upper()} lidar: NEES/NIS inside 99% bounds and MSE at the CRB",
              f1 > 0.90 and f2 > 0.90 and 0.85 < eff < 1.2,
              f"(NEES {f1:.1%}, NIS {f2:.1%}, MSE/CRB {eff:.3f})")


def test_burn_makes_range_data_limited():
    """Coast: range bound scales with the prior. After a burn it does not.
    Same truth noise (same seeds) for both cases, averaged over 10 truths."""
    out = {}
    for burn in (False, True):
        kw = {"burn_time": 0.5 * 5578.2, "burn_dv": (0.1, 0, 0)} if burn else {}
        sc = Scenario(lidar=False, **kw)
        weak = Scenario(lidar=False, **kw)
        weak.P0 = 25.0 * sc.P0
        r1, r5 = [], []
        for seed in range(10):
            X, _ = sc.truth(np.random.default_rng(seed))
            r1.append(range_sigma(X[-1], info_bound(sc, X)[-1]))
            r5.append(range_sigma(X[-1], info_bound(weak, X)[-1]))
        out[burn] = (np.mean(r1), np.mean(r5))
    c, b = out[False], out[True]
    check("angles-only coast: range bound grows with the prior (5x sigma)", c[1] / c[0] > 4.0,
          f"({c[0]:.0f} -> {c[1]:.0f} m)")
    check("after a 0.1 m/s radial burn the range bound barely depends on the prior",
          b[1] / b[0] < 1.3 and b[0] < 0.6 * c[0], f"({b[0]:.0f} -> {b[1]:.0f} m)")


if __name__ == "__main__":
    test_stm()
    test_van_loan()
    test_jacobian()
    test_angles_only_observability()
    test_ukf_equals_kf_linear()
    test_dare()
    test_rts_equals_batch()
    test_burn_makes_range_data_limited()
    test_monte_carlo_consistency()
    print("all checks passed")
