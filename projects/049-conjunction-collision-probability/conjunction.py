# Astron Built this AP
# Project 049 | Orbital mechanics | Conjunction collision probability and avoidance burn sizing
# Author: Aditya (aditya-astron)
# Date: 2026-11-20 IST
"""
Short-term encounter collision probability for two objects at TCA.

For a high relative velocity the encounter lasts milliseconds, so relative
motion is rectilinear and the positional error is frozen. Projecting onto the
plane normal to the relative velocity (the encounter or B-plane) reduces the
problem to a 2D Gaussian integrated over a disk of the combined hard-body
radius R:

    Pc = integral_{|x| <= R} N(x; m, C) dA

with m the projected miss vector and C the projected combined covariance.
Four evaluations are provided, each with a different failure mode:

  pc_foster  polar quadrature over the disk (Gauss-Legendre in r, trapezoid in
             theta, which is spectrally accurate for a periodic integrand)
  pc_1d      exact reduction to a 1D integral of normal CDFs in the principal
             frame (Alfano), integrated adaptively
  pc_chan    Chan's series with the equivalent-area isotropic substitution;
             exact for an isotropic C, approximate otherwise
  pc_small   small-object limit, pi R^2 times the density at the disk centre
"""

import math

import numpy as np
from scipy import integrate
from scipy.special import ndtr

MU_EARTH = 3.986004418e14  # m^3/s^2
R_EARTH = 6378137.0        # m


# ---------------------------------------------------------------------------
# Frames
# ---------------------------------------------------------------------------

def rtn_basis(r, v):
    """Columns are the radial, transverse, normal unit vectors in ECI."""
    r_hat = r / np.linalg.norm(r)
    h = np.cross(r, v)
    n_hat = h / np.linalg.norm(h)
    t_hat = np.cross(n_hat, r_hat)
    return np.column_stack([r_hat, t_hat, n_hat])


def rtn_to_eci_cov(r, v, cov_rtn):
    q = rtn_basis(r, v)
    return q @ cov_rtn @ q.T


def encounter_frame(dr, dv):
    """Rows are the encounter-frame axes in ECI.

    z along the relative velocity, x along the miss vector (the component of
    dr normal to dv), y completes the right-handed set. The first two rows
    span the encounter plane.
    """
    z = dv / np.linalg.norm(dv)
    x = dr - (dr @ z) * z
    x = x / np.linalg.norm(x)
    y = np.cross(z, x)
    return np.vstack([x, y, z])


def project_encounter(dr, dv, cov3):
    """Projected 2D miss vector and covariance in the encounter plane."""
    e = encounter_frame(dr, dv)
    m = (e @ dr)[:2]
    c = (e @ cov3 @ e.T)[:2, :2]
    return m, 0.5 * (c + c.T)


def rectilinear_miss(dr, dv):
    """Time of closest approach and miss distance for straight-line motion."""
    t = -(dr @ dv) / (dv @ dv)
    return t, np.linalg.norm(dr + dv * t)


# ---------------------------------------------------------------------------
# Probability of collision
# ---------------------------------------------------------------------------

def _principal(m, c):
    """Rotate (m, C) to principal axes: returns sigma_x >= sigma_y and m'."""
    w, v = np.linalg.eigh(c)
    w, v = w[::-1], v[:, ::-1]
    if np.any(w <= 0.0):
        raise ValueError("projected covariance is not positive definite")
    return np.sqrt(w), v.T @ m


def gauss2(x, y, m, c):
    """Bivariate normal density, vectorised over x and y."""
    ci = np.linalg.inv(c)
    dx, dy = x - m[0], y - m[1]
    q = ci[0, 0] * dx * dx + 2.0 * ci[0, 1] * dx * dy + ci[1, 1] * dy * dy
    return np.exp(-0.5 * q) / (2.0 * math.pi * math.sqrt(np.linalg.det(c)))


def pc_foster(m, c, radius, n_r=48, n_theta=96):
    """Polar quadrature over the hard-body disk."""
    xg, wg = np.polynomial.legendre.leggauss(n_r)
    r = 0.5 * radius * (xg + 1.0)
    wr = 0.5 * radius * wg
    th = 2.0 * math.pi * np.arange(n_theta) / n_theta
    rr, tt = np.meshgrid(r, th, indexing="ij")
    f = gauss2(rr * np.cos(tt), rr * np.sin(tt), m, c)
    return float((wr * r) @ f.sum(axis=1) * (2.0 * math.pi / n_theta))


def pc_1d(m, c, radius):
    """Alfano's exact 1D form, x = R sin(t) removes the chord-end singularity.

    Pc = int_{-R}^{R} phi((x-mx)/sx)/sx * [Phi((h-my)/sy) - Phi((-h-my)/sy)] dx,
    h = sqrt(R^2 - x^2), in the principal frame of C.
    """
    (sx, sy), (mx, my) = _principal(m, c)

    def f(t):
        x = radius * math.sin(t)
        h = radius * math.cos(t)
        px = math.exp(-0.5 * ((x - mx) / sx) ** 2) / (sx * math.sqrt(2.0 * math.pi))
        hi, lo = (h - my) / sy, (-h - my) / sy
        # Difference of CDFs taken in the tail it lives in, so a chord far
        # out on the positive side does not cancel to 1 - 1.
        band = ndtr(hi) - ndtr(lo) if hi + lo < 0.0 else ndtr(-lo) - ndtr(-hi)
        return px * band * h

    # Concentrate effort where the density peaks inside the chord range.
    t_pk = math.asin(max(-1.0, min(1.0, mx / radius)))
    val, _ = integrate.quad(f, -0.5 * math.pi, 0.5 * math.pi, points=[t_pk],
                            epsabs=0.0, epsrel=1e-13, limit=400)
    return val


def pc_chan(m, c, radius, tol=1e-16, max_terms=2000):
    """Chan's series with u = R^2/(sx sy), v = mx^2/sx^2 + my^2/sy^2.

    Pc = exp(-v/2) sum_m (v/2)^m/m! [1 - exp(-u/2) sum_{k<=m} (u/2)^k/k!]

    which is the noncentral chi-square CDF P(chi'^2_2(v) <= u). Exact when
    C is isotropic; otherwise the ellipse is replaced by a circle of equal
    Mahalanobis area.
    """
    (sx, sy), (mx, my) = _principal(m, c)
    u = radius * radius / (sx * sy)
    v = (mx / sx) ** 2 + (my / sy) ** 2
    a, b = 0.5 * u, 0.5 * v
    pois_v = math.exp(-b)   # e^{-v/2} (v/2)^m / m!
    pois_u = math.exp(-a)   # e^{-u/2} (u/2)^k / k!
    cdf_u = pois_u          # running sum_{k<=m}
    total, m_term = 0.0, 0
    while True:
        total += pois_v * (1.0 - cdf_u)
        m_term += 1
        pois_v *= b / m_term
        pois_u *= a / m_term
        cdf_u += pois_u
        if (pois_v * (1.0 - cdf_u) < tol * max(total, 1e-300) and m_term > b) \
                or m_term >= max_terms:
            return total


def pc_small(m, c, radius):
    """Small-object limit: density at the disk centre times its area."""
    return math.pi * radius * radius * float(gauss2(0.0, 0.0, m, c))


# ---------------------------------------------------------------------------
# Covariance scaling (dilution)
# ---------------------------------------------------------------------------

def pc_scaled(m, c, radius, k):
    """Pc when both standard deviations are scaled by k (C -> k^2 C)."""
    return pc_foster(m, k * k * c, radius)


def max_pc_small_object(m, c, radius):
    """Peak of the small-object Pc over C -> k^2 C.

    Pc(k) = R^2/(2 k^2 sqrt|C|) exp(-d^2/(2 k^2)), d^2 = m' C^-1 m, has its
    maximum at k^2 = d^2/2 with Pc_max = R^2 / (e d^2 sqrt|C|).
    """
    d2 = float(m @ np.linalg.solve(c, m))
    return math.sqrt(0.5 * d2), radius ** 2 / (math.e * d2 * math.sqrt(np.linalg.det(c)))


def max_pc_numeric(m, c, radius, k_lo=1e-2, k_hi=1e2):
    """Golden-section search on log k for the covariance-scaling maximum."""
    from scipy.optimize import minimize_scalar
    res = minimize_scalar(lambda lk: -pc_scaled(m, c, radius, math.exp(lk)),
                          bounds=(math.log(k_lo), math.log(k_hi)),
                          method="bounded", options={"xatol": 1e-10})
    return math.exp(res.x), -res.fun


# ---------------------------------------------------------------------------
# Monte Carlo on the 3D rectilinear encounter
# ---------------------------------------------------------------------------

def pc_monte_carlo(dr0, dv0, cov_r, cov_v, radius, n, rng, chunk=500_000):
    """Sample relative position and velocity, count straight-line misses < R.

    Independent of the encounter-plane projection: it samples in ECI and
    finds each sample's own closest approach and direction of motion.
    """
    lr = np.linalg.cholesky(cov_r)
    lv = np.linalg.cholesky(cov_v)
    hits, done = 0, 0
    while done < n:
        k = min(chunk, n - done)
        dr = dr0 + rng.standard_normal((k, 3)) @ lr.T
        dv = dv0 + rng.standard_normal((k, 3)) @ lv.T
        t = -np.einsum("ij,ij->i", dr, dv) / np.einsum("ij,ij->i", dv, dv)
        miss = np.linalg.norm(dr + dv * t[:, None], axis=1)
        hits += int(np.count_nonzero(miss < radius))
        done += k
    return hits / n, hits


# ---------------------------------------------------------------------------
# Avoidance manoeuvre: CW prediction and two-body truth
# ---------------------------------------------------------------------------

def cw_impulse_displacement(n, t, dv_rtn):
    """Clohessy-Wiltshire position at time t after an impulse dv (R, T, N)."""
    vx, vy, vz = dv_rtn
    s, c = math.sin(n * t), math.cos(n * t)
    x = s / n * vx + 2.0 / n * (1.0 - c) * vy
    y = -2.0 / n * (1.0 - c) * vx + (4.0 * s / n - 3.0 * t) * vy
    z = s / n * vz
    return np.array([x, y, z])


def _twobody(t, y):
    r = y[:3]
    return np.concatenate([y[3:], -MU_EARTH * r / np.linalg.norm(r) ** 3])


def propagate(state, dt):
    if dt == 0.0:
        return state.copy()
    sol = integrate.solve_ivp(_twobody, (0.0, dt), state, method="DOP853",
                              rtol=1e-13, atol=1e-7)
    return sol.y[:, -1]


def twobody_impulse_displacement(r_tca, v_tca, lead, dv_rtn):
    """Displacement at TCA (in TCA RTN) of an impulse applied `lead` s before.

    Both nominal and manoeuvred arcs start from the same back-propagated
    state, so integration error cancels in the difference.
    """
    s0 = propagate(np.concatenate([r_tca, v_tca]), -lead)
    q0 = rtn_basis(s0[:3], s0[3:])
    s1 = s0.copy()
    s1[3:] += q0 @ np.asarray(dv_rtn, dtype=float)
    a = propagate(s0, lead)
    b = propagate(s1, lead)
    return rtn_basis(r_tca, v_tca).T @ (b[:3] - a[:3])
