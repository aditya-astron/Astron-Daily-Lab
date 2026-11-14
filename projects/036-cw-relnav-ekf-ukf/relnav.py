# Astron Built this AP
# Project 036 | Kalman filtering and sensor fusion | Rendezvous relative navigation in the Hill frame: EKF vs UKF vs the PCRB
# Author: Aditya (aditya-astron)
# Date: 2026-11-14 IST
"""
Relative navigation for a chaser closing on a target in a near-circular LEO.

State x = [x, y, z, vx, vy, vz] of the chaser in the target's Hill (LVLH)
frame: x radial out, y along-track, z orbit normal. Clohessy-Wiltshire
dynamics with white unmodelled acceleration (differential drag, thruster
leakage) of PSD q per axis:

    x'' = 3 n^2 x + 2 n y' + w_x
    y'' =        - 2 n x' + w_y
    z'' =   -n^2 z        + w_z

The chaser sees the target along u = -r. A scanning lidar returns
[range, azimuth, elevation]; a camera returns only the two angles

    az = atan2(u_x, u_y),   el = asin(u_z / |u|).

Two cases:

  A. Lidar on a coelliptic drift from 8 km behind. EKF and UKF are checked
     for statistical consistency (Monte Carlo NEES/NIS against chi-square
     bounds) and for efficiency against the posterior Cramer-Rao bound.

  B. Camera only. Angles are homogeneous of degree zero in r, and the CW
     dynamics are linear, so a trajectory scaled by any alpha > 0 produces
     exactly the same angles: range is unobservable on a coast. A known
     impulse breaks the scaling (Woffinden & Geller 2009). The PCRB and a
     UKF Monte Carlo show the range collapse after a small radial burn.
"""

import numpy as np
from scipy.linalg import cholesky

MU_EARTH = 3.986004418e14       # m^3/s^2
R_EARTH = 6378137.0             # m


def mean_motion(alt_m):
    a = R_EARTH + alt_m
    return np.sqrt(MU_EARTH / a**3)


# --------------------------------------------------------------------------
# Dynamics
# --------------------------------------------------------------------------

def cw_A(n):
    A = np.zeros((6, 6))
    A[0:3, 3:6] = np.eye(3)
    A[3, 0] = 3 * n * n
    A[3, 4] = 2 * n
    A[4, 3] = -2 * n
    A[5, 2] = -n * n
    return A


def cw_phi(n, t):
    """Closed-form CW state transition matrix."""
    s, c, nt = np.sin(n * t), np.cos(n * t), n * t
    P = np.zeros((6, 6))
    P[0, 0] = 4 - 3 * c
    P[1, 0] = 6 * (s - nt)
    P[1, 1] = 1.0
    P[2, 2] = c
    P[0, 3] = s / n
    P[0, 4] = 2 * (1 - c) / n
    P[1, 3] = -2 * (1 - c) / n
    P[1, 4] = (4 * s - 3 * nt) / n
    P[2, 5] = s / n
    P[3, 0] = 3 * n * s
    P[4, 0] = -6 * n * (1 - c)
    P[5, 2] = -n * s
    P[3, 3] = c
    P[3, 4] = 2 * s
    P[4, 3] = -2 * s
    P[4, 4] = 4 * c - 3
    P[5, 5] = c
    return P


def van_loan_Q(n, dt, q):
    """Discrete process noise Q = int_0^dt Phi(s) G q G^T Phi(s)^T ds.

    Van Loan: expm of [[-A, GqG^T], [0, A^T]] dt gives Phi^-1 Q in the
    upper-right block and Phi^T in the lower-right.
    """
    from scipy.linalg import expm
    A = cw_A(n)
    G = np.vstack([np.zeros((3, 3)), np.eye(3)])
    M = np.zeros((12, 12))
    M[:6, :6] = -A
    M[:6, 6:] = G @ (q * np.eye(3)) @ G.T
    M[6:, 6:] = A.T
    E = expm(M * dt)
    PhiT = E[6:, 6:]
    Q = PhiT.T @ E[:6, 6:]
    return 0.5 * (Q + Q.T)


# --------------------------------------------------------------------------
# Measurements
# --------------------------------------------------------------------------

def h_lidar(x):
    u = -x[:3]
    rho = np.linalg.norm(u)
    return np.array([rho, np.arctan2(u[0], u[1]), np.arcsin(u[2] / rho)])


def H_lidar(x):
    u = -x[:3]
    ux, uy, uz = u
    rho2 = u @ u
    rho = np.sqrt(rho2)
    rxy2 = ux * ux + uy * uy
    rxy = np.sqrt(rxy2)
    J = np.zeros((3, 3))            # d[rho, az, el] / du
    J[0] = u / rho
    J[1] = [uy / rxy2, -ux / rxy2, 0.0]
    J[2] = [-ux * uz / (rho2 * rxy), -uy * uz / (rho2 * rxy), rxy / rho2]
    H = np.zeros((3, 6))
    H[:, :3] = -J                   # du/dr = -I
    return H


def h_angles(x):
    return h_lidar(x)[1:]


def H_angles(x):
    return H_lidar(x)[1:]


def wrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def innovation(z, zhat, angle_idx):
    nu = z - zhat
    nu[angle_idx] = wrap(nu[angle_idx])
    return nu


# --------------------------------------------------------------------------
# Filters. Dynamics are linear, so both share the exact prediction; they
# differ only in how the nonlinear measurement update is approximated.
# --------------------------------------------------------------------------

def predict(x, P, Phi, Q):
    return Phi @ x, Phi @ P @ Phi.T + Q


def ekf_update(x, P, z, h, Hf, R, angle_idx):
    H = Hf(x)
    nu = innovation(z, h(x), angle_idx)
    S = H @ P @ H.T + R
    K = np.linalg.solve(S, H @ P).T
    IKH = np.eye(len(x)) - K @ H
    P = IKH @ P @ IKH.T + K @ R @ K.T          # Joseph form
    return x + K @ nu, 0.5 * (P + P.T), nu, S


class UT:
    """Scaled unscented transform weights (Wan & van der Merwe)."""

    def __init__(self, n, alpha=1.0, beta=2.0, kappa=0.0):
        lam = alpha**2 * (n + kappa) - n
        self.c = n + lam
        self.wm = np.full(2 * n + 1, 0.5 / self.c)
        self.wc = self.wm.copy()
        self.wm[0] = lam / self.c
        self.wc[0] = lam / self.c + 1 - alpha**2 + beta

    def sigma(self, x, P):
        L = cholesky(self.c * P, lower=True)
        return np.vstack([x, x + L.T, x - L.T])


def ukf_update(x, P, z, h, R, angle_idx, ut):
    X = ut.sigma(x, P)
    Z = np.array([h(s) for s in X])
    # Angle means: average residuals about the central point, not raw angles.
    dZ = Z - Z[0]
    dZ[:, angle_idx] = wrap(dZ[:, angle_idx])
    zhat = Z[0] + ut.wm @ dZ
    zhat[angle_idx] = wrap(zhat[angle_idx])
    dZ = Z - zhat
    dZ[:, angle_idx] = wrap(dZ[:, angle_idx])
    dX = X - x
    S = (ut.wc * dZ.T) @ dZ + R
    Pxz = (ut.wc * dX.T) @ dZ
    K = np.linalg.solve(S, Pxz.T).T
    nu = innovation(z, zhat, angle_idx)
    P = P - K @ S @ K.T
    return x + K @ nu, 0.5 * (P + P.T), nu, S


def rts_smooth(xf, Pf, xp, Pp, Phi):
    """Rauch-Tung-Striebel smoother over stored filter/predict sequences."""
    N = len(xf)
    xs, Ps = xf.copy(), Pf.copy()
    for k in range(N - 2, -1, -1):
        C = np.linalg.solve(Pp[k + 1], Phi @ Pf[k]).T
        xs[k] = xf[k] + C @ (xs[k + 1] - xp[k + 1])
        Ps[k] = Pf[k] + C @ (Ps[k + 1] - Pp[k + 1]) @ C.T
    return xs, Ps


# --------------------------------------------------------------------------
# Scenario
# --------------------------------------------------------------------------

class Scenario:
    def __init__(self, alt=420e3, dt=10.0, orbits=1.5, q=1e-10,
                 sig_rho=1.0, sig_ang=1e-3, lidar=True,
                 burn_time=None, burn_dv=(0.0, 0.0, 0.0)):
        self.n = mean_motion(alt)
        self.T = 2 * np.pi / self.n
        self.dt = dt
        self.N = int(round(orbits * self.T / dt))
        self.Phi = cw_phi(self.n, dt)
        self.Q = van_loan_Q(self.n, dt, q)
        self.Lq = cholesky(self.Q, lower=True)
        self.lidar = lidar
        if lidar:
            self.h, self.Hf = h_lidar, H_lidar
            self.R = np.diag([sig_rho**2, sig_ang**2, sig_ang**2])
            self.ang = [1, 2]
        else:
            self.h, self.Hf = h_angles, H_angles
            self.R = np.diag([sig_ang**2, sig_ang**2])
            self.ang = [0, 1]
        self.Lr = np.sqrt(np.diag(self.R))
        # Coelliptic drift 500 m below and 8 km behind: vy = 1.5 n |dx|.
        self.x0 = np.array([-500.0, -8000.0, 100.0,
                            0.0, 1.5 * self.n * 500.0, 0.0])
        # Hand-over from relative GPS / ground tracking.
        self.P0 = np.diag([100.0**2, 500.0**2, 50.0**2,
                           0.02**2, 0.05**2, 0.01**2])
        # Burn applied at the end of step k_burn (known to the filter).
        self.k_burn = None if burn_time is None else int(round(burn_time / dt))
        self.dv = np.r_[0.0, 0.0, 0.0, np.asarray(burn_dv, float)]

    def propagate_mean(self, x, k):
        x = self.Phi @ x
        if self.k_burn is not None and k == self.k_burn:
            x = x + self.dv
        return x

    def truth(self, rng):
        X = np.empty((self.N + 1, 6))
        Z = np.empty((self.N + 1, len(self.Lr)))
        X[0] = self.x0 + cholesky(self.P0, lower=True) @ rng.standard_normal(6)
        for k in range(self.N):
            X[k + 1] = self.propagate_mean(X[k], k) + self.Lq @ rng.standard_normal(6)
        for k in range(self.N + 1):
            Z[k] = self.h(X[k]) + self.Lr * rng.standard_normal(len(self.Lr))
            Z[k, self.ang] = wrap(Z[k, self.ang])
        return X, Z

    def run_filter(self, Z, kind, ut=None, store_pred=False):
        """Return filtered states, covariances and NIS at every step."""
        N = self.N
        xf = np.empty((N + 1, 6))
        Pf = np.empty((N + 1, 6, 6))
        nis = np.empty(N + 1)
        if store_pred:
            xp = np.empty((N + 1, 6))
            Pp = np.empty((N + 1, 6, 6))
        x, P = self.x0.copy(), self.P0.copy()
        if ut is None and kind == "ukf":
            ut = UT(6)
        for k in range(N + 1):
            if k > 0:
                x, P = predict(x, P, self.Phi, self.Q)
                if self.k_burn is not None and k - 1 == self.k_burn:
                    x = x + self.dv
            if store_pred:
                xp[k], Pp[k] = x, P
            if kind == "ekf":
                x, P, nu, S = ekf_update(x, P, Z[k], self.h, self.Hf, self.R, self.ang)
            else:
                x, P, nu, S = ukf_update(x, P, Z[k], self.h, self.R, self.ang, ut)
            nis[k] = nu @ np.linalg.solve(S, nu)
            xf[k], Pf[k] = x, P
        if store_pred:
            return xf, Pf, nis, xp, Pp
        return xf, Pf, nis


def info_bound(sc, X):
    """Cramer-Rao bound along one true trajectory X (prior P0, linear
    dynamics, additive Gaussian noise):

        J_{k+1} = (Q + Phi J_k^-1 Phi^T)^-1 + H(x_{k+1})^T R^-1 H(x_{k+1})

    This is the information filter run on the truth instead of the estimate.
    Returns the bound covariance J^-1 at every step.
    """
    Ri = np.linalg.inv(sc.R)
    J = np.linalg.inv(sc.P0)
    out = np.empty((sc.N + 1, 6, 6))
    for k in range(sc.N + 1):
        if k > 0:
            J = np.linalg.inv(sc.Q + sc.Phi @ np.linalg.solve(J, sc.Phi.T))
        H = sc.Hf(X[k])
        J = J + H.T @ Ri @ H
        out[k] = np.linalg.inv(J)
    return out


def pcrb_bayes(sc, truths):
    """Bayesian posterior CRB (Tichavsky et al. 1998): same recursion with
    E[H^T R^-1 H] averaged over the Monte Carlo truths. It is a valid lower
    bound on the MSE but loose by Jensen's inequality when truths spread in
    range, and it mixes the null spaces of different trajectories, so it
    hides the angles-only range unobservability that info_bound exposes.
    """
    Ri = np.linalg.inv(sc.R)
    J = np.linalg.inv(sc.P0)
    out = np.empty((sc.N + 1, 6, 6))
    for k in range(sc.N + 1):
        if k > 0:
            J = np.linalg.inv(sc.Q + sc.Phi @ np.linalg.solve(J, sc.Phi.T))
        J = J + sum(sc.Hf(X[k]).T @ Ri @ sc.Hf(X[k]) for X in truths) / len(truths)
        out[k] = np.linalg.inv(J)
    return out


def monte_carlo(sc, runs, seed, kinds=("ekf", "ukf")):
    rng = np.random.default_rng(seed)
    truths, res = [], {k: {"err": [], "P": [], "nis": []} for k in kinds}
    for _ in range(runs):
        X, Z = sc.truth(rng)
        truths.append(X)
        for kind in kinds:
            xf, Pf, nis = sc.run_filter(Z, kind)
            res[kind]["err"].append(xf - X)
            res[kind]["P"].append(Pf)
            res[kind]["nis"].append(nis)
    for kind in kinds:
        e = np.array(res[kind]["err"])
        P = np.array(res[kind]["P"])
        res[kind]["err"] = e
        res[kind]["P"] = P
        res[kind]["nis"] = np.array(res[kind]["nis"])
        res[kind]["nees"] = np.einsum("rki,rki->rk", e, np.linalg.solve(P, e[..., None])[..., 0])
    return truths, res


def range_sigma(x, P):
    """1-sigma of range |r| from the position covariance (first order)."""
    r = x[:3]
    g = r / np.linalg.norm(r)
    return np.sqrt(g @ P[:3, :3] @ g)


def mean_bound(sc, truths):
    """Average over truths of the per-trajectory bound covariance."""
    return np.mean([info_bound(sc, X) for X in truths], axis=0)


def range_error(truths, err):
    Xt = np.array(truths)
    return (np.linalg.norm(Xt[..., :3] + err[..., :3], axis=-1)
            - np.linalg.norm(Xt[..., :3], axis=-1))


def main():
    from scipy.stats import chi2

    print("Rendezvous relative navigation, target at 420 km circular")
    runs = 100

    # ---------------- Case A: lidar ----------------
    sc = Scenario(lidar=True)
    print(f"  n = {sc.n:.6e} rad/s, T = {sc.T:.1f} s, dt = {sc.dt:.0f} s, "
          f"{sc.N} steps ({sc.N * sc.dt / sc.T:.2f} orbits)")
    xe = sc.x0.copy()
    for k in range(sc.N):
        xe = sc.Phi @ xe
    print(f"  nominal coelliptic drift: r0 = {sc.x0[:3]} m -> r_end = {np.round(xe[:3], 1)} m")
    print("\nCase A: lidar (1 m, 1 mrad), q = 1e-10 m^2/s^3, "
          f"{runs} Monte Carlo runs")
    truths, res = monte_carlo(sc, runs, seed=1)
    Bt = mean_bound(sc, truths)
    Bb = pcrb_bayes(sc, truths)
    lo, hi = chi2.ppf([0.005, 0.995], runs * 6) / runs
    lo3, hi3 = chi2.ppf([0.005, 0.995], runs * 3) / runs
    tail = slice(sc.N // 3, None)
    print(f"  99% bounds on run-averaged NEES [{lo:.2f}, {hi:.2f}], NIS [{lo3:.2f}, {hi3:.2f}]")
    print("  filter  NEES in  NIS in  mean NEES  mean NIS  pos RMSE end  vel RMSE end  MSE/bound")
    bt = np.trace(Bt[:, :3, :3], axis1=1, axis2=2)
    for kind in ("ekf", "ukf"):
        r = res[kind]
        anees, anis, e = r["nees"].mean(0), r["nis"].mean(0), r["err"]
        mse_pos = (e[:, :, :3] ** 2).sum(-1).mean(0)
        print(f"  {kind.upper():6s} {np.mean((anees > lo) & (anees < hi)):7.1%} "
              f"{np.mean((anis > lo3) & (anis < hi3)):7.1%} {anees.mean():10.2f} "
              f"{anis.mean():9.2f} {np.sqrt(mse_pos[-1]):11.2f} m "
              f"{1e3 * np.sqrt((e[:, -1, 3:] ** 2).sum(-1).mean()):8.2f} mm/s "
              f"{np.mean(mse_pos[tail] / bt[tail]):9.3f}")
    print(f"  per-trajectory CRB at end: pos {np.sqrt(bt[-1]):.2f} m, "
          f"vel {1e3 * np.sqrt(np.trace(Bt[-1, 3:, 3:])):.2f} mm/s "
          f"(MSE/bound averaged over the last 2/3 of the arc)")
    print(f"  Bayesian PCRB at end: pos {np.sqrt(np.trace(Bb[-1, :3, :3])):.2f} m "
          "(looser: averages information over the range spread)")

    rng = np.random.default_rng(7)
    X, Z = sc.truth(rng)
    xf, Pf, _, xp, Pp = sc.run_filter(Z, "ekf", store_pred=True)
    xs, Ps = rts_smooth(xf, Pf, xp, Pp, sc.Phi)
    print(f"  RTS smoother, one run: hand-over state pos sigma "
          f"{np.sqrt(np.trace(Pf[0, :3, :3])):.1f} m filtered -> "
          f"{np.sqrt(np.trace(Ps[0, :3, :3])):.2f} m smoothed "
          f"(actual error {np.linalg.norm(xs[0, :3] - X[0, :3]):.2f} m)")

    # ---------------- Case B: angles only ----------------
    print(f"\nCase B: camera angles only (1 mrad), {runs} Monte Carlo runs")
    tb = 0.5 * sc.T
    dv = (0.10, 0.0, 0.0)
    for label, kw in (("coast", {}),
                      (f"radial burn {dv[0]:.2f} m/s at t = {tb:.0f} s",
                       {"burn_time": tb, "burn_dv": dv})):
        sb = Scenario(lidar=False, **kw)
        truths_b, res_b = monte_carlo(sb, runs, seed=2)
        bounds = [info_bound(sb, X) for X in truths_b]
        weak = Scenario(lidar=False, **kw)
        weak.P0 = 25.0 * sb.P0
        bounds_w = [info_bound(weak, X) for X in truths_b]
        print(f"  {label}:")
        for name, k in (("t = T/2 + 100 s", sb.N // 2 + 10), ("end", sb.N)):
            rs = np.mean([range_sigma(X[k], B[k]) for X, B in zip(truths_b, bounds)])
            rw = np.mean([range_sigma(X[k], B[k]) for X, B in zip(truths_b, bounds_w)])
            msg = f"    {name:16s} CRB range sigma {rs:6.1f} m (5x prior: {rw:6.1f} m)"
            for kind in ("ekf", "ukf"):
                re = range_error(truths_b, res_b[kind]["err"])[:, k]
                msg += f" | {kind.upper()} RMSE {np.sqrt(np.mean(re ** 2)):6.1f} m"
            print(msg)
        print("    mean NEES (ideal 6): " + ", ".join(
            f"{kind.upper()} {res_b[kind]['nees'].mean():.2f}" for kind in ("ekf", "ukf")))


if __name__ == "__main__":
    main()
