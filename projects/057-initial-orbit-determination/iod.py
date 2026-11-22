# Astron Built this AP
# Project 057 | Orbital mechanics / Navigation | Initial orbit determination: Gibbs, Herrick-Gibbs and Gauss angles-only
# Author: Aditya (aditya-astron)
# Date: 2026-11-22 IST
"""Initial orbit determination (IOD) from three sightings.

Batch least squares (project 026) needs a starting state. IOD supplies it
with no prior: three position vectors (radar) give the velocity by Gibbs or
Herrick-Gibbs, three line-of-sight directions (telescope) give the whole
state by Gauss's method. Everything here is two-body; the exact propagator
is the universal-variable Kepler solution, which also feeds the Gauss
refinement step.

Units: km, s, rad. Formulation follows Curtis (Orbital Mechanics for
Engineering Students, ch. 3 and 5) and Vallado (Fundamentals of
Astrodynamics and Applications, ch. 7).
"""
from __future__ import annotations

import numpy as np

MU = 398600.4418          # km^3/s^2
RE = 6378.137             # km, WGS84 equatorial radius
F_EARTH = 1.0 / 298.257223563
OMEGA_E = 7.292115e-5     # rad/s


# ---------------------------------------------------------------------------
# Universal-variable two-body propagation
# ---------------------------------------------------------------------------

def stumpff_c(z: float) -> float:
    # Series near z = 0 avoids the 1 - cos cancellation that loses ~8 digits.
    if abs(z) < 1e-3:
        return 0.5 - z / 24.0 + z * z / 720.0 - z ** 3 / 40320.0
    if z > 0:
        return (1.0 - np.cos(np.sqrt(z))) / z
    return (np.cosh(np.sqrt(-z)) - 1.0) / (-z)


def stumpff_s(z: float) -> float:
    if abs(z) < 1e-3:
        return 1.0 / 6.0 - z / 120.0 + z * z / 5040.0 - z ** 3 / 362880.0
    if z > 0:
        sz = np.sqrt(z)
        return (sz - np.sin(sz)) / sz ** 3
    sz = np.sqrt(-z)
    return (np.sinh(sz) - sz) / sz ** 3


def universal_anomaly(dt: float, r0: float, vr0: float, alpha: float,
                      mu: float = MU, tol: float = 1e-12, maxit: int = 60) -> float:
    """Newton solve of the universal Kepler equation for chi at time dt.

    F(chi) = r0 vr0/sqrt(mu) chi^2 C(z) + (1 - alpha r0) chi^3 S(z)
             + r0 chi - sqrt(mu) dt,          z = alpha chi^2
    dF/dchi equals the radius r(chi), so the Newton step is (F / r).
    """
    smu = np.sqrt(mu)
    chi = smu * abs(alpha) * dt if alpha > 1e-12 else smu * dt / r0
    for _ in range(maxit):
        z = alpha * chi * chi
        C, S = stumpff_c(z), stumpff_s(z)
        F = (r0 * vr0 / smu * chi * chi * C + (1.0 - alpha * r0) * chi ** 3 * S
             + r0 * chi - smu * dt)
        dF = (r0 * vr0 / smu * chi * (1.0 - z * S)
              + (1.0 - alpha * r0) * chi * chi * C + r0)
        step = F / dF
        chi -= step
        if abs(step) < tol * max(1.0, abs(chi)):
            return chi
    raise RuntimeError("universal Kepler iteration did not converge")


def lagrange_fg(r0v: np.ndarray, v0v: np.ndarray, dt: float, mu: float = MU):
    """Exact Lagrange coefficients f, g, fdot, gdot over dt."""
    r0 = np.linalg.norm(r0v)
    vr0 = r0v @ v0v / r0
    alpha = 2.0 / r0 - v0v @ v0v / mu
    chi = universal_anomaly(dt, r0, vr0, alpha, mu)
    z = alpha * chi * chi
    C, S = stumpff_c(z), stumpff_s(z)
    f = 1.0 - chi * chi / r0 * C
    g = dt - chi ** 3 / np.sqrt(mu) * S
    r = np.linalg.norm(f * r0v + g * v0v)
    fdot = np.sqrt(mu) / (r * r0) * (alpha * chi ** 3 * S - chi)
    gdot = 1.0 - chi * chi / r * C
    return f, g, fdot, gdot


def propagate(r0v, v0v, dt, mu: float = MU):
    r0v, v0v = np.asarray(r0v, float), np.asarray(v0v, float)
    f, g, fdot, gdot = lagrange_fg(r0v, v0v, dt, mu)
    return f * r0v + g * v0v, fdot * r0v + gdot * v0v


# ---------------------------------------------------------------------------
# Elements
# ---------------------------------------------------------------------------

def coe_to_rv(a, e, i, raan, argp, nu, mu: float = MU):
    p = a * (1 - e * e)
    r_pf = p / (1 + e * np.cos(nu)) * np.array([np.cos(nu), np.sin(nu), 0.0])
    v_pf = np.sqrt(mu / p) * np.array([-np.sin(nu), e + np.cos(nu), 0.0])
    cO, sO, ci, si, cw, sw = (np.cos(raan), np.sin(raan), np.cos(i), np.sin(i),
                              np.cos(argp), np.sin(argp))
    Q = np.array([
        [cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci, sO * si],
        [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si],
        [sw * si, cw * si, ci]])
    return Q @ r_pf, Q @ v_pf


def rv_to_coe(rv, vv, mu: float = MU):
    """Return (a, e, i, raan, argp, nu) for an elliptic, inclined orbit."""
    r, v = np.linalg.norm(rv), np.linalg.norm(vv)
    h = np.cross(rv, vv)
    n = np.cross([0.0, 0.0, 1.0], h)
    ev = ((v * v - mu / r) * rv - (rv @ vv) * vv) / mu
    e = np.linalg.norm(ev)
    a = 1.0 / (2.0 / r - v * v / mu)
    i = np.arccos(h[2] / np.linalg.norm(h))
    raan = np.arctan2(n[1], n[0]) % (2 * np.pi)
    argp = np.arccos(np.clip(n @ ev / (np.linalg.norm(n) * e), -1, 1))
    if ev[2] < 0:
        argp = 2 * np.pi - argp
    nu = np.arccos(np.clip(ev @ rv / (e * r), -1, 1))
    if rv @ vv < 0:
        nu = 2 * np.pi - nu
    return a, e, i, raan, argp, nu


# ---------------------------------------------------------------------------
# Velocity from three position vectors
# ---------------------------------------------------------------------------

def coplanarity(r1, r2, r3) -> float:
    """Angle (rad) of r1 out of the plane of r2, r3. Three two-body
    positions are coplanar; a large value means bad data, not an orbit."""
    c23 = np.cross(r2, r3)
    return np.arcsin(abs(r1 @ c23) / (np.linalg.norm(r1) * np.linalg.norm(c23)))


def gibbs(r1, r2, r3, mu: float = MU) -> np.ndarray:
    """Gibbs' method: purely geometric, no times needed.

    N = r1 (r2 x r3) + r2 (r3 x r1) + r3 (r1 x r2)
    D = r2 x r3 + r3 x r1 + r1 x r2
    S = r1 (r2 - r3) + r2 (r3 - r1) + r3 (r1 - r2)
    v2 = sqrt(mu / (|N| |D|)) (D x r2 / |r2| + S)

    All three are differences of nearly equal vectors when the arc is
    short, which is why Gibbs degrades on closely spaced data.
    """
    n1, n2, n3 = (np.linalg.norm(x) for x in (r1, r2, r3))
    c12, c23, c31 = np.cross(r1, r2), np.cross(r2, r3), np.cross(r3, r1)
    N = n1 * c23 + n2 * c31 + n3 * c12
    D = c12 + c23 + c31
    S = r1 * (n2 - n3) + r2 * (n3 - n1) + r3 * (n1 - n2)
    return np.sqrt(mu / (np.linalg.norm(N) * np.linalg.norm(D))) * (
        np.cross(D, r2) / n2 + S)


def herrick_gibbs(r1, r2, r3, t1, t2, t3, mu: float = MU) -> np.ndarray:
    """Herrick-Gibbs: Taylor-series velocity with the two-body acceleration
    folded in. Truncation error is O(dt^4) in velocity, so it wins over
    Gibbs on short arcs and loses on long ones.
    """
    d21, d31, d32 = t2 - t1, t3 - t1, t3 - t2
    n1, n2, n3 = (np.linalg.norm(x) for x in (r1, r2, r3))
    return (-d32 * (1.0 / (d21 * d31) + mu / (12.0 * n1 ** 3)) * r1
            + (d32 - d21) * (1.0 / (d21 * d32) + mu / (12.0 * n2 ** 3)) * r2
            + d21 * (1.0 / (d32 * d31) + mu / (12.0 * n3 ** 3)) * r3)


# ---------------------------------------------------------------------------
# Observation geometry
# ---------------------------------------------------------------------------

def site_eci(lat: float, alt: float, theta: float,
             re: float = RE, flat: float = F_EARTH) -> np.ndarray:
    """Inertial position of a site on the oblate Earth (Curtis eq. 5.56);
    theta is the local sidereal time (Earth rotation angle + east longitude)."""
    e2 = 2 * flat - flat * flat
    den = np.sqrt(1.0 - e2 * np.sin(lat) ** 2)
    rxy = (re / den + alt) * np.cos(lat)
    rz = (re * (1 - flat) ** 2 / den + alt) * np.sin(lat)
    return np.array([rxy * np.cos(theta), rxy * np.sin(theta), rz])


def radec_to_unit(ra: float, dec: float) -> np.ndarray:
    return np.array([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def elevation(r_obj, r_site) -> float:
    """Elevation above the local horizontal of an oblate site, using the
    geodetic normal recovered from the site vector's latitude."""
    rho = r_obj - r_site
    # Up direction on the ellipsoid: gradient of x^2+y^2 + z^2/(1-f)^2.
    up = np.array([r_site[0], r_site[1], r_site[2] / (1 - F_EARTH) ** 2])
    up /= np.linalg.norm(up)
    return np.arcsin(rho @ up / np.linalg.norm(rho))


# ---------------------------------------------------------------------------
# Gauss angles-only IOD
# ---------------------------------------------------------------------------

def gauss(rho_hat, R, t, mu: float = MU, refine: str | None = "fixed-point",
          tol: float = 1e-9, maxit: int = 200, root: str = "auto"):
    """Gauss's method for three topocentric unit vectors rho_hat[k] seen from
    inertial site positions R[k] at times t[k].

    Returns dict(r2, v2, rho, r2_initial, v2_initial, roots, iterations).

    refine: None keeps the series answer, "fixed-point" is Curtis's
    iteration, "correct" hands the series answer to angles_correct().

    The middle range comes from the 8th-degree polynomial
        x^8 + a x^6 + b x^3 + c = 0     (x = |r2|)
    built from the series f, g truncated after the mu/r^3 term. The
    refinement then replaces the series coefficients with exact universal
    f, g and re-solves the linear range equations until the ranges stop
    moving (Curtis algorithm 5.6).
    """
    L1, L2, L3 = (np.asarray(x, float) for x in rho_hat)
    R1, R2, R3 = (np.asarray(x, float) for x in R)
    tau1, tau3 = t[0] - t[1], t[2] - t[1]
    tau = tau3 - tau1

    p1, p2, p3 = np.cross(L2, L3), np.cross(L1, L3), np.cross(L1, L2)
    D0 = L1 @ p1
    if abs(D0) < 1e-14:
        raise ValueError("lines of sight are coplanar (D0 = 0); Gauss is singular")
    Dm = np.array([[Ri @ pj for pj in (p1, p2, p3)] for Ri in (R1, R2, R3)])
    # Dm[i, j] = R_(i+1) . p_(j+1), matching Curtis's D_ij.

    A = (-Dm[0, 1] * tau3 / tau + Dm[1, 1] + Dm[2, 1] * tau1 / tau) / D0
    B = (Dm[0, 1] * (tau3 ** 2 - tau ** 2) * tau3 / tau
         + Dm[2, 1] * (tau ** 2 - tau1 ** 2) * tau1 / tau) / (6.0 * D0)
    E = R2 @ L2
    a = -(A * A + 2 * A * E + R2 @ R2)
    b = -2 * mu * B * (A + E)
    c = -(mu * B) ** 2
    coeffs = [1, 0, a, 0, 0, b, 0, 0, c]
    rts = np.roots(coeffs)
    pos = sorted(float(z.real) for z in rts if abs(z.imag) < 1e-6 * abs(z) and z.real > 0)
    if not pos:
        raise ValueError("no positive real root of the Gauss polynomial")
    if root == "auto":
        # Reject roots that put the object below the surface or behind the
        # observer; of what remains take the one with the smallest positive
        # range (the usual single-solution case leaves exactly one).
        cand = [x for x in pos if x > RE and A + mu * B / x ** 3 > 0]
        if not cand:
            raise ValueError("no physical root of the Gauss polynomial")
        r2m = cand[0]
    else:
        r2m = pos[int(root)]

    r3c = r2m ** 3
    rho2 = A + mu * B / r3c
    rho1 = ((6 * (Dm[2, 0] * tau1 / tau3 + Dm[1, 0] * tau / tau3) * r3c
             + mu * Dm[2, 0] * (tau ** 2 - tau1 ** 2) * tau1 / tau3)
            / (6 * r3c + mu * (tau ** 2 - tau3 ** 2)) - Dm[0, 0]) / D0
    rho3 = ((6 * (Dm[0, 2] * tau3 / tau1 - Dm[1, 2] * tau / tau1) * r3c
             + mu * Dm[0, 2] * (tau ** 2 - tau3 ** 2) * tau3 / tau1)
            / (6 * r3c + mu * (tau ** 2 - tau1 ** 2)) - Dm[2, 2]) / D0

    r1, r2, r3 = R1 + rho1 * L1, R2 + rho2 * L2, R3 + rho3 * L3
    f1 = 1 - 0.5 * mu * tau1 ** 2 / r3c
    f3 = 1 - 0.5 * mu * tau3 ** 2 / r3c
    g1 = tau1 - mu * tau1 ** 3 / (6 * r3c)
    g3 = tau3 - mu * tau3 ** 3 / (6 * r3c)
    v2 = (-f3 * r1 + f1 * r3) / (f1 * g3 - f3 * g1)
    out = dict(r2_initial=r2.copy(), v2_initial=v2.copy(), roots=pos, iterations=0)

    if refine == "correct":
        r2, v2, it = angles_correct(rho_hat, R, t, r2, v2, mu)
        out.update(r2=r2, v2=v2, iterations=it,
                   rho=np.array([np.linalg.norm(propagate(r2, v2, tk - t[1], mu)[0] - Rk)
                                 for tk, Rk in zip(t, (R1, R2, R3))]))
        return out
    if refine == "fixed-point":
        rho_old = np.array([rho1, rho2, rho3])
        for it in range(1, maxit + 1):
            f1n, g1n, _, _ = lagrange_fg(r2, v2, tau1, mu)
            f3n, g3n, _, _ = lagrange_fg(r2, v2, tau3, mu)
            # Averaging with the previous coefficients damps the fixed-point
            # iteration; it converges on long arcs where the plain update
            # oscillates.
            f1, f3, g1, g3 = (f1 + f1n) / 2, (f3 + f3n) / 2, (g1 + g1n) / 2, (g3 + g3n) / 2
            den = f1 * g3 - f3 * g1
            c1, c3 = g3 / den, -g1 / den
            rho1 = (-Dm[0, 0] + Dm[1, 0] / c1 - Dm[2, 0] * c3 / c1) / D0
            rho2 = (-c1 * Dm[0, 1] + Dm[1, 1] - c3 * Dm[2, 1]) / D0
            rho3 = (-Dm[0, 2] * c1 / c3 + Dm[1, 2] / c3 - Dm[2, 2]) / D0
            r1, r2, r3 = R1 + rho1 * L1, R2 + rho2 * L2, R3 + rho3 * L3
            v2 = (-f3 * r1 + f1 * r3) / den
            rho_new = np.array([rho1, rho2, rho3])
            if np.max(np.abs(rho_new - rho_old)) < tol * np.max(np.abs(rho_new)):
                out["iterations"] = it
                break
            rho_old = rho_new
        else:
            raise RuntimeError("Gauss refinement did not converge")

    out.update(r2=r2, v2=v2, rho=np.array([rho1, rho2, rho3]))
    return out


def _los_residual(L, R, t, x, mu):
    """Six residuals: each predicted sightline projected on two unit vectors
    normal to the observed one. Zero exactly when all three match."""
    res = np.empty(6)
    for k in range(3):
        r = propagate(x[:3], x[3:], t[k] - t[1], mu)[0] if k != 1 else x[:3]
        u = (r - R[k]) / np.linalg.norm(r - R[k])
        e1 = np.cross(L[k], [0.0, 0.0, 1.0])
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(L[k], e1)
        res[2 * k:2 * k + 2] = (u @ e1, u @ e2)
    return res


def angles_correct(rho_hat, R, t, r2, v2, mu: float = MU,
                   tol: float = 1e-10, maxit: int = 50):
    """Newton differential correction of (r2, v2) on the three sightlines.

    Three angle pairs fix the six-element state exactly, so this converges
    to the same orbit as a converged Gauss fixed point; unlike that
    iteration it keeps converging on long arcs (tens of degrees), where
    the f, g update overshoots. Jacobian by central differences, step
    halving if a full step increases the residual.
    """
    L = [np.asarray(x, float) / np.linalg.norm(x) for x in rho_hat]
    R = [np.asarray(x, float) for x in R]
    x = np.concatenate([r2, v2]).astype(float)
    h = np.array([1e-3] * 3 + [1e-6] * 3)
    res = _los_residual(L, R, t, x, mu)
    for it in range(1, maxit + 1):
        J = np.empty((6, 6))
        for j in range(6):
            dx = np.zeros(6)
            dx[j] = h[j]
            J[:, j] = (_los_residual(L, R, t, x + dx, mu)
                       - _los_residual(L, R, t, x - dx, mu)) / (2 * h[j])
        step = np.linalg.solve(J, -res)
        lam = 1.0
        while True:
            xn = x + lam * step
            try:
                rn = _los_residual(L, R, t, xn, mu)
                ok = np.linalg.norm(rn) < np.linalg.norm(res) or lam < 1e-3
            except RuntimeError:          # step left the elliptic regime badly
                ok = False
            if ok:
                break
            lam /= 2
        x, res = xn, rn
        if np.linalg.norm(lam * step[:3]) < tol * np.linalg.norm(x[:3]):
            return x[:3], x[3:], it
    raise RuntimeError("angles-only differential correction did not converge")
