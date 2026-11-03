# Astron Built this AP
# Project 009 | Kalman filtering and sensor fusion | Multiplicative EKF fusing a drifting gyro with a star tracker
# Author: Aditya (aditya-astron)
# Date: 2026-11-03 IST
"""Multiplicative extended Kalman filter (MEKF) for spacecraft attitude.

A rate-integrating gyro is propagated at 10 Hz and corrected by a star tracker
quaternion at 1 Hz. The filter carries a global unit quaternion plus a 6-state
error vector x = [dtheta, dbeta]: a small body-frame rotation error and the
gyro bias error. Keeping the attitude error three-dimensional avoids the
singular 4x4 covariance an additive quaternion EKF would have.

Conventions
  q = [w, x, y, z], Hamilton product, rotates body vectors into inertial.
  Kinematics  q_dot = 0.5 q (x) [0, w_body]
  Error       q_true = q_est (x) dq,  dq ~ [1, dtheta/2]

Gyro model (Farrenkopf):
  w_meas = w_true + beta + eta_v        eta_v: angle random walk, sigma_v
  beta_dot = eta_u                      eta_u: rate random walk,  sigma_u

Error dynamics with w_hat = w_meas - beta_hat:
  dtheta_dot = -[w_hat x] dtheta - dbeta - eta_v
  dbeta_dot  = eta_u
"""
import numpy as np
from scipy.linalg import expm, solve_discrete_are

ARCSEC = np.pi / 180.0 / 3600.0
DEG_PER_HR = np.pi / 180.0 / 3600.0


# --- quaternion helpers -----------------------------------------------------

def qmul(p, q):
    pw, pv = p[0], p[1:]
    qw, qv = q[0], q[1:]
    return np.concatenate(([pw * qw - pv @ qv], pw * qv + qw * pv + np.cross(pv, qv)))


def qconj(q):
    return np.array([q[0], -q[1], -q[2], -q[3]])


def qexp(phi):
    """Quaternion for a rotation vector phi (rad)."""
    a = np.linalg.norm(phi)
    if a < 1e-12:
        q = np.concatenate(([1.0], 0.5 * phi))
        return q / np.linalg.norm(q)
    return np.concatenate(([np.cos(a / 2)], np.sin(a / 2) * phi / a))


def qerr(q_est, q_true):
    """Body-frame rotation vector taking q_est onto q_true."""
    dq = qmul(qconj(q_est), q_true)
    if dq[0] < 0:
        dq = -dq
    s = np.linalg.norm(dq[1:])
    if s < 1e-15:
        return 2.0 * dq[1:]
    return 2.0 * np.arctan2(s, dq[0]) * dq[1:] / s


def skew(v):
    return np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])


# --- discretisation ---------------------------------------------------------

def discretize(w_hat, dt, sigma_v, sigma_u):
    """Exact Phi and Qd for the error model via Van Loan's matrix exponential."""
    F = np.zeros((6, 6))
    F[:3, :3] = -skew(w_hat)
    F[:3, 3:] = -np.eye(3)
    G = np.zeros((6, 6))
    G[:3, :3] = -np.eye(3)
    G[3:, 3:] = np.eye(3)
    Qc = np.diag([sigma_v**2] * 3 + [sigma_u**2] * 3)

    M = np.zeros((12, 12))
    M[:6, :6] = -F
    M[:6, 6:] = G @ Qc @ G.T
    M[6:, 6:] = F.T
    E = expm(M * dt)
    Phi = E[6:, 6:].T
    Qd = Phi @ E[:6, 6:]
    return Phi, 0.5 * (Qd + Qd.T)


def farrenkopf_qd(dt, sigma_v, sigma_u):
    """Closed-form Qd for zero body rate (Crassidis & Junkins, Optimal Estimation)."""
    I = np.eye(3)
    q11 = (sigma_v**2 * dt + sigma_u**2 * dt**3 / 3.0) * I
    q12 = -(sigma_u**2 * dt**2 / 2.0) * I
    q22 = sigma_u**2 * dt * I
    return np.block([[q11, q12], [q12, q22]])


# --- filter -----------------------------------------------------------------

class MEKF:
    def __init__(self, q0, beta0, P0, sigma_v, sigma_u, sigma_st):
        self.q = q0 / np.linalg.norm(q0)
        self.beta = np.array(beta0, dtype=float)
        self.P = np.array(P0, dtype=float)
        self.sigma_v = sigma_v
        self.sigma_u = sigma_u
        self.R = sigma_st**2 * np.eye(3)
        self.H = np.hstack([np.eye(3), np.zeros((3, 3))])

    def propagate(self, w_meas, dt):
        w_hat = w_meas - self.beta
        self.q = qmul(self.q, qexp(w_hat * dt))
        self.q /= np.linalg.norm(self.q)
        Phi, Qd = discretize(w_hat, dt, self.sigma_v, self.sigma_u)
        self.P = Phi @ self.P @ Phi.T + Qd

    def update(self, q_meas):
        # Innovation is the small rotation from the estimate to the measurement.
        y = qerr(self.q, q_meas)
        S = self.H @ self.P @ self.H.T + self.R
        K = np.linalg.solve(S, self.H @ self.P).T
        dx = K @ y
        # Reset: fold the attitude correction into the global quaternion.
        self.q = qmul(self.q, qexp(dx[:3]))
        self.q /= np.linalg.norm(self.q)
        self.beta = self.beta + dx[3:]
        # Joseph form keeps P symmetric positive definite over long runs.
        IKH = np.eye(6) - K @ self.H
        self.P = IKH @ self.P @ IKH.T + K @ self.R @ K.T
        return y, S


# --- truth and sensor simulation -------------------------------------------

def truth_rate(t):
    """Slow tumble with a sinusoidal slew superimposed, rad/s."""
    d = np.pi / 180.0
    return np.array([0.10 * d * np.sin(0.01 * t),
                     0.05 * d * np.cos(0.007 * t),
                     0.08 * d])


def simulate(T=3600.0, dt=0.1, st_period=1.0, sigma_v=np.sqrt(10) * 1e-7,
             sigma_u=np.sqrt(10) * 1e-10, sigma_st=6.0 * ARCSEC,
             beta_true0=None, seed=0, rate_fn=truth_rate):
    """Run truth, sensors and filter. Returns a dict of time histories."""
    rng = np.random.default_rng(seed)
    if beta_true0 is None:
        beta_true0 = np.array([0.1, -0.1, 0.15]) * DEG_PER_HR
    n = int(round(T / dt))
    st_every = int(round(st_period / dt))

    q_true = qexp(np.array([0.3, -0.2, 0.5]))
    beta = np.array(beta_true0, dtype=float)

    # Filter starts 0.1 deg off with zero bias knowledge.
    P0 = np.diag([(0.1 * np.pi / 180) ** 2] * 3 + [(0.2 * DEG_PER_HR) ** 2] * 3)
    q0 = qmul(q_true, qexp(rng.normal(0, 0.1 * np.pi / 180, 3)))
    f = MEKF(q0, np.zeros(3), P0, sigma_v, sigma_u, sigma_st)

    t_hist, err, sig, nees = [], [], [], []
    for k in range(1, n + 1):
        t = k * dt
        # Integrate truth at the midpoint rate (piecewise constant over dt).
        w = rate_fn(t - 0.5 * dt)
        beta_prev = beta.copy()
        beta = beta + sigma_u * np.sqrt(dt) * rng.standard_normal(3)
        # Rate-integrating gyro: mean bias over the interval plus ARW and the
        # integrated RRW term, so the measured angle has Farrenkopf statistics.
        w_meas = (w + 0.5 * (beta + beta_prev)
                  + np.sqrt(sigma_v**2 / dt + sigma_u**2 * dt / 12.0) * rng.standard_normal(3))
        q_true = qmul(q_true, qexp(w * dt))
        q_true /= np.linalg.norm(q_true)

        f.propagate(w_meas, dt)
        if k % st_every == 0:
            q_meas = qmul(q_true, qexp(sigma_st * rng.standard_normal(3)))
            f.update(q_meas)
            e = np.concatenate((qerr(f.q, q_true), beta - f.beta))
            t_hist.append(t)
            err.append(e)
            sig.append(np.sqrt(np.diag(f.P)))
            nees.append(e @ np.linalg.solve(f.P, e))

    return {"t": np.array(t_hist), "err": np.array(err), "sig": np.array(sig),
            "nees": np.array(nees), "P": f.P.copy()}


# --- analytic steady state -------------------------------------------------

def single_axis_steady_state(dt_meas, sigma_v, sigma_u, sigma_st):
    """Steady-state posterior [attitude var, bias var] for one axis at zero rate.

    Between star tracker updates the gyro propagation is a pure integration, so
    the many small gyro steps combine into one exact step of length dt_meas.
    Solve the discrete algebraic Riccati equation for the prior and apply one
    measurement update.
    """
    A = np.array([[1.0, -dt_meas], [0.0, 1.0]])
    Q = farrenkopf_qd(dt_meas, sigma_v, sigma_u)[[0, 3]][:, [0, 3]]
    H = np.array([[1.0, 0.0]])
    R = np.array([[sigma_st**2]])
    Pm = solve_discrete_are(A.T, H.T, Q, R)
    K = Pm @ H.T / (H @ Pm @ H.T + R)
    Pp = (np.eye(2) - K @ H) @ Pm
    return Pp[0, 0], Pp[1, 1]


def main():
    sv, su, sn = np.sqrt(10) * 1e-7, np.sqrt(10) * 1e-10, 6.0 * ARCSEC
    r = simulate(T=3600.0, seed=1, sigma_v=sv, sigma_u=su, sigma_st=sn)
    ss = r["t"] > 1200.0
    rms_att = np.sqrt(np.mean(r["err"][ss, :3] ** 2, axis=0)) / ARCSEC
    rms_b = np.sqrt(np.mean(r["err"][ss, 3:] ** 2, axis=0)) / DEG_PER_HR
    p_att, p_b = single_axis_steady_state(1.0, sv, su, sn)

    print("MEKF: 10 Hz gyro + 1 Hz star tracker, 1 h, slewing at ~0.14 deg/s")
    print(f"  gyro ARW {sv:.3e} rad/s^0.5, RRW {su:.3e} rad/s^1.5, ST {sn / ARCSEC:.1f} arcsec")
    print("  steady state (t > 1200 s)       x        y        z    analytic")
    print("  attitude RMS [arcsec]     " + " ".join(f"{v:8.3f}" for v in rms_att)
          + f"  {np.sqrt(p_att) / ARCSEC:8.3f}")
    print("  filter 1-sigma [arcsec]   " + " ".join(f"{v:8.3f}" for v in r["sig"][-1, :3] / ARCSEC))
    print("  bias RMS [deg/hr]         " + " ".join(f"{v:8.4f}" for v in rms_b)
          + f"  {np.sqrt(p_b) / DEG_PER_HR:8.4f}")
    print("  filter 1-sigma [deg/hr]   " + " ".join(f"{v:8.4f}" for v in r["sig"][-1, 3:] / DEG_PER_HR))
    print(f"  mean NEES (6 dof)         {r['nees'][ss].mean():.2f}")
    print(f"  attitude improvement over raw star tracker: {sn / np.sqrt(p_att):.2f}x")


if __name__ == "__main__":
    main()
