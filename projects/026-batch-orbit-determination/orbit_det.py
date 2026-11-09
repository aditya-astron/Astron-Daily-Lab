# Astron Built this AP
# Project 026 | Navigation | Batch least-squares orbit determination from ground-station range and range-rate
# Author: Aditya (aditya-astron)
# Date: 2026-11-09 IST
"""
Batch weighted least-squares orbit determination of a LEO spacecraft.

Three ground stations (Bengaluru, Hartebeesthoek, Wallops) record two-way range and
range-rate whenever the spacecraft is above a 10 deg elevation mask. From a
24 h arc we estimate the epoch state x0 = [r0, v0] (ECI, km and km/s).

Dynamics: two-body + J2,

    a = -mu r / |r|^3 + (3/2) J2 mu Re^2 / |r|^5 * [x(5s-1), y(5s-1), z(5s-3)],
    s = z^2 / |r|^2

State transition matrix Phi(t, t0) = d x(t) / d x0 is integrated alongside the
state from  dPhi/dt = A(t) Phi,  A = [[0, I], [G, 0]],  G = da/dr (analytic).

Measurements, with station position R and velocity Vs = w_E x R in ECI:

    rho    = |d|,           d = r - R
    rhodot = d . (v - Vs) / rho

Partials (H~ = dG/dx(t)):

    d rho / d r     = d / rho                       d rho / d v    = 0
    d rhodot / d r  = w / rho - (d.w) d / rho^3     d rhodot / d v = d / rho

with w = v - Vs. The batch maps them to the epoch with H = H~ Phi(t, t0) and
solves the Gauss-Newton normal equations

    (H^T W H) dx0 = H^T W dy,   W = diag(1/sigma^2),   P0 = (H^T W H)^-1

iterating until the update is negligible. Light time, troposphere, station
clock and Earth orientation beyond a uniform spin are not modelled; the same
model generates the truth, so the filter's formal covariance must match the
scatter seen in Monte Carlo.
"""
import numpy as np
from scipy.integrate import solve_ivp

MU = 398600.4418          # km^3/s^2
RE = 6378.137             # km, WGS-84 equatorial radius
J2 = 1.08262668e-3
W_E = 7.2921150e-5        # rad/s
F_WGS = 1 / 298.257223563
THETA0 = 0.0              # Earth rotation angle at epoch, rad

SIG_RHO = 10e-3           # km    (10 m)
SIG_RDOT = 1e-6           # km/s  (1 mm/s)
EL_MASK = np.radians(10.0)

STATIONS = {              # name: (geodetic lat deg, lon deg, height km)
    "Bengaluru": (13.034, 77.512, 0.85),
    "Hartebeesthoek": (-25.887, 27.707, 1.56),
    "Wallops": (37.940, -75.466, 0.01),
}


# ---------------------------------------------------------------- dynamics
def accel(r):
    x, y, z = r
    rn2 = r @ r
    rn = np.sqrt(rn2)
    s = z * z / rn2
    k = 1.5 * J2 * MU * RE**2 / rn**5
    return -MU * r / rn**3 + k * np.array([x * (5 * s - 1), y * (5 * s - 1), z * (5 * s - 3)])


def gravity_gradient(r):
    """G = da/dr, 3x3, two-body + J2 (analytic)."""
    x, y, z = r
    rn2 = r @ r
    rn = np.sqrt(rn2)
    G = -MU / rn**3 * (np.eye(3) - 3 * np.outer(r, r) / rn2)

    # J2 term: a = k' f, k' = 1.5 J2 mu Re^2,
    # f_x = x g, f_y = y g, f_z = z h, g = 5 z^2 r^-7 - r^-5, h = g - 2 r^-5
    kp = 1.5 * J2 * MU * RE**2
    r5, r7, r9 = rn**-5, rn**-7, rn**-9
    g = 5 * z * z * r7 - r5
    h = 5 * z * z * r7 - 3 * r5
    ez = np.array([0.0, 0.0, 1.0])
    dg = 10 * z * r7 * ez + (-35 * z * z * r9 + 5 * r7) * r
    dh = 10 * z * r7 * ez + (-35 * z * z * r9 + 15 * r7) * r
    GJ = np.empty((3, 3))
    GJ[0] = x * dg + np.array([g, 0, 0])
    GJ[1] = y * dg + np.array([0, g, 0])
    GJ[2] = z * dh + np.array([0, 0, h])
    return G + kp * GJ


def _rhs(t, Y):
    r, v = Y[:3], Y[3:6]
    Phi = Y[6:].reshape(6, 6)
    A = np.zeros((6, 6))
    A[:3, 3:] = np.eye(3)
    A[3:, :3] = gravity_gradient(r)
    return np.concatenate([v, accel(r), (A @ Phi).ravel()])


def propagate(x0, times, rtol=1e-11, atol=1e-12):
    """State and STM at each time in `times` (seconds from epoch, sorted, times[0] >= 0)."""
    Y0 = np.concatenate([x0, np.eye(6).ravel()])
    t_eval = np.asarray(times, float)
    sol = solve_ivp(_rhs, (0.0, t_eval[-1]), Y0, method="DOP853",
                    t_eval=t_eval, rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    X = sol.y[:6].T
    Phi = sol.y[6:].T.reshape(-1, 6, 6)
    return X, Phi


# ---------------------------------------------------------------- stations
def geodetic_to_ecef(lat_deg, lon_deg, h_km):
    lat, lon = np.radians(lat_deg), np.radians(lon_deg)
    e2 = F_WGS * (2 - F_WGS)
    N = RE / np.sqrt(1 - e2 * np.sin(lat) ** 2)
    return np.array([(N + h_km) * np.cos(lat) * np.cos(lon),
                     (N + h_km) * np.cos(lat) * np.sin(lon),
                     (N * (1 - e2) + h_km) * np.sin(lat)])


def station_eci(R_ecef, t):
    th = THETA0 + W_E * t
    c, s = np.cos(th), np.sin(th)
    R = np.array([c * R_ecef[0] - s * R_ecef[1], s * R_ecef[0] + c * R_ecef[1], R_ecef[2]])
    V = W_E * np.array([-R[1], R[0], 0.0])
    return R, V


def up_vector(lat_deg, lon_deg, t):
    """Geodetic local vertical in ECI."""
    lat, lon = np.radians(lat_deg), np.radians(lon_deg) + THETA0 + W_E * t
    return np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])


# ---------------------------------------------------------------- measurements
def range_rangerate(x, R, V):
    d = x[:3] - R
    w = x[3:] - V
    rho = np.linalg.norm(d)
    return rho, d @ w / rho


def meas_partials(x, R, V):
    """2x6 Jacobian of [rho, rhodot] w.r.t. the state at measurement time."""
    d = x[:3] - R
    w = x[3:] - V
    rho = np.linalg.norm(d)
    dw = d @ w
    Ht = np.zeros((2, 6))
    Ht[0, :3] = d / rho
    Ht[1, :3] = w / rho - dw * d / rho**3
    Ht[1, 3:] = d / rho
    return Ht


def schedule(x_true, t_end, dt):
    """Visible (time, station) pairs over the arc under the elevation mask."""
    times = np.arange(0.0, t_end + 0.5 * dt, dt)
    X, _ = propagate(x_true, times)
    obs = []
    for i, t in enumerate(times):
        for name, (lat, lon, h) in STATIONS.items():
            R, _ = station_eci(geodetic_to_ecef(lat, lon, h), t)
            d = X[i, :3] - R
            el = np.arcsin(d @ up_vector(lat, lon, t) / np.linalg.norm(d))
            if el > EL_MASK:
                obs.append((t, name))
    return obs


def simulate(x_true, obs, rng=None):
    """Noisy (or exact, rng=None) [rho, rhodot] for each scheduled observation."""
    times = np.array(sorted({t for t, _ in obs}))
    X, _ = propagate(x_true, times)
    idx = {t: i for i, t in enumerate(times)}
    Y = np.empty((len(obs), 2))
    for j, (t, name) in enumerate(obs):
        R, V = station_eci(geodetic_to_ecef(*STATIONS[name]), t)
        Y[j] = range_rangerate(X[idx[t]], R, V)
    if rng is not None:
        Y += rng.standard_normal(Y.shape) * [SIG_RHO, SIG_RDOT]
    return Y


# ---------------------------------------------------------------- estimator
def batch_ls(x_guess, obs, Y, max_iter=20, tol=1e-6, verbose=False):
    """Gauss-Newton batch estimate of the epoch state. Returns x0, P0, info."""
    times = np.array(sorted({t for t, _ in obs}))
    idx = {t: i for i, t in enumerate(times)}
    stn = {n: geodetic_to_ecef(*STATIONS[n]) for n in STATIONS}
    w = np.array([1 / SIG_RHO**2, 1 / SIG_RDOT**2])
    # Scale velocity columns so the normal matrix is well conditioned (km vs km/s).
    S = np.diag([1, 1, 1, 1e3, 1e3, 1e3]).astype(float)
    x = np.array(x_guess, float)
    history = []
    for it in range(1, max_iter + 1):
        Xs, Phis = propagate(x, times)
        N = np.zeros((6, 6))
        b = np.zeros(6)
        res = np.empty_like(Y)
        for j, (t, name) in enumerate(obs):
            R, V = station_eci(stn[name], t)
            xi = Xs[idx[t]]
            res[j] = Y[j] - np.array(range_rangerate(xi, R, V))
            H = meas_partials(xi, R, V) @ Phis[idx[t]] @ S
            N += H.T @ (w[:, None] * H)
            b += H.T @ (w * res[j])
        dz = np.linalg.solve(N, b)
        dx = S @ dz
        cost = float(np.sum(res**2 * w))
        history.append((it, cost, np.linalg.norm(dx[:3]), np.linalg.norm(dx[3:])))
        if verbose:
            print(f"  iter {it}: cost {cost:14.4f}  |dr| {np.linalg.norm(dx[:3])*1e3:12.4f} m"
                  f"  |dv| {np.linalg.norm(dx[3:])*1e6:12.4f} mm/s")
        x = x + dx
        if np.linalg.norm(dx[:3]) < tol and np.linalg.norm(dx[3:]) < tol * 1e-3:
            break
    else:
        raise RuntimeError("batch did not converge")
    # Final residuals and covariance at the converged state.
    Xs, Phis = propagate(x, times)
    N = np.zeros((6, 6))
    res = np.empty_like(Y)
    for j, (t, name) in enumerate(obs):
        R, V = station_eci(stn[name], t)
        res[j] = Y[j] - np.array(range_rangerate(Xs[idx[t]], R, V))
        H = meas_partials(Xs[idx[t]], R, V) @ Phis[idx[t]] @ S
        N += H.T @ (w[:, None] * H)
    P = S @ np.linalg.inv(N) @ S
    info = {"iterations": len(history), "history": history, "residuals": res,
            "cost": float(np.sum(res**2 * w)), "dof": Y.size - 6}
    return x, P, info


def rtn_frame(x):
    r, v = x[:3], x[3:]
    R = r / np.linalg.norm(r)
    N = np.cross(r, v)
    N /= np.linalg.norm(N)
    return np.vstack([R, np.cross(N, R), N])


def coe_to_rv(a, e, i, raan, argp, nu):
    p = a * (1 - e * e)
    rp = p / (1 + e * np.cos(nu)) * np.array([np.cos(nu), np.sin(nu), 0])
    vp = np.sqrt(MU / p) * np.array([-np.sin(nu), e + np.cos(nu), 0])
    cO, sO, ci, si, cw, sw = (np.cos(raan), np.sin(raan), np.cos(i), np.sin(i),
                              np.cos(argp), np.sin(argp))
    Q = np.array([[cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci, sO * si],
                  [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si],
                  [sw * si, cw * si, ci]])
    return np.concatenate([Q @ rp, Q @ vp])


# 500 km x 51.6 deg, slightly eccentric
X_TRUE = coe_to_rv(RE + 500.0, 0.001, np.radians(51.6), np.radians(40.0),
                   np.radians(30.0), np.radians(10.0))
ARC_S = 24 * 3600.0
DT_S = 10.0


def main():
    rng = np.random.default_rng(26)
    obs = schedule(X_TRUE, ARC_S, DT_S)
    per = {n: sum(1 for _, s in obs if s == n) for n in STATIONS}
    print(f"Arc {ARC_S/3600:.0f} h, {DT_S:.0f} s sampling, {len(obs)} range/range-rate pairs")
    for n, c in per.items():
        print(f"  {n:<15s} {c:4d} pairs")
    Y = simulate(X_TRUE, obs, rng)

    x0 = X_TRUE + np.array([2.0, -1.5, 1.0, 2e-3, -1e-3, 1.5e-3])   # ~2.7 km, 2.7 m/s off
    print("\nGauss-Newton iterations:")
    x_hat, P, info = batch_ls(x0, obs, Y, verbose=True)

    res = info["residuals"]
    print(f"\nPost-fit residual RMS: range {np.sqrt(np.mean(res[:,0]**2))*1e3:.2f} m "
          f"(sigma {SIG_RHO*1e3:.0f} m), range-rate {np.sqrt(np.mean(res[:,1]**2))*1e6:.3f} mm/s "
          f"(sigma {SIG_RDOT*1e6:.0f} mm/s)")
    print(f"Weighted cost {info['cost']:.1f} for {info['dof']} dof "
          f"(reduced chi-square {info['cost']/info['dof']:.3f})")

    T = rtn_frame(X_TRUE)
    err = x_hat - X_TRUE
    sr = np.sqrt(np.diag(T @ P[:3, :3] @ T.T))
    sv = np.sqrt(np.diag(T @ P[3:, 3:] @ T.T))
    er, ev = T @ err[:3], T @ err[3:]
    print("\nEpoch state error vs formal 1-sigma (RTN):")
    for k, lab in enumerate("RTN"):
        print(f"  {lab}  pos {er[k]*1e3:8.3f} m  (1s {sr[k]*1e3:6.3f} m)   "
              f"vel {ev[k]*1e6:8.4f} mm/s  (1s {sv[k]*1e6:6.4f} mm/s)")
    nees = err @ np.linalg.solve(P, err)
    print(f"  NEES (6 dof, expect ~6): {nees:.2f}")


if __name__ == "__main__":
    main()
