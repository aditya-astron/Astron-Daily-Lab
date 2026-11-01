# Astron Built this AP
# Project 003 | Orbital mechanics | Lambert targeting for a LEO rendezvous transfer
# Author: Aditya (aditya-astron)
# Date: 2026-11-01 IST
"""Single-revolution Lambert solver (universal variables) and a two-body
propagator used to check its answers independently."""

import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

MU_EARTH = 398600.4418  # km^3/s^2
R_EARTH = 6378.137      # km


# Closed forms cancel badly near z = 0 (1 - cos(sqrt z) at z = 1e-8 keeps
# only ~8 digits), so switch to the Taylor series for |z| < 0.1.
#   C(z) = sum (-z)^k / (2k+2)!,   S(z) = sum (-z)^k / (2k+3)!
# Six terms leave a remainder below 1e-17 at |z| = 0.1.
_C_COEF = [1 / math.factorial(2 * k + 2) for k in range(6)]
_S_COEF = [1 / math.factorial(2 * k + 3) for k in range(6)]
_Z_SERIES = 0.1


def _series(coef, z):
    acc = 0.0
    for c in reversed(coef):
        acc = acc * (-z) + c
    return acc


def stumpff_c(z):
    if z > _Z_SERIES:
        s = np.sqrt(z)
        return (1.0 - np.cos(s)) / z
    if z < -_Z_SERIES:
        s = np.sqrt(-z)
        return (np.cosh(s) - 1.0) / (-z)
    return _series(_C_COEF, z)


def stumpff_s(z):
    if z > _Z_SERIES:
        s = np.sqrt(z)
        return (s - np.sin(s)) / s**3
    if z < -_Z_SERIES:
        s = np.sqrt(-z)
        return (np.sinh(s) - s) / s**3
    return _series(_S_COEF, z)


def lambert(r1, r2, tof, mu=MU_EARTH, prograde=True, min_sin=2e-3):
    """Velocities (v1, v2) of the conic joining r1 -> r2 in time tof.

    Universal-variable form (Bate-Mueller-White / Curtis Alg. 5.2):
        A    = sin(dth) * sqrt(r1 r2 / (1 - cos dth))
        y(z) = r1 + r2 + A (z S(z) - 1) / sqrt(C(z))
        F(z) = (y / C)^1.5 S + A sqrt(y) - sqrt(mu) tof = 0
    z = alpha * chi^2 sets the conic: z < 0 hyperbola, z > 0 ellipse.
    F is monotonic in z on the single-rev branch, so a bracketed root finder
    is safer than the usual Newton iteration.
    """
    r1 = np.asarray(r1, float)
    r2 = np.asarray(r2, float)
    n1, n2 = np.linalg.norm(r1), np.linalg.norm(r2)

    cos_dth = np.clip(np.dot(r1, r2) / (n1 * n2), -1.0, 1.0)
    dth = np.arccos(cos_dth)
    hz = np.cross(r1, r2)[2]
    if (prograde and hz < 0) or (not prograde and hz >= 0):
        dth = 2.0 * np.pi - dth
    # At 0/180 deg the plane is undefined and g -> 0; the f-g velocities lose
    # digits fast (LEO case: propagated miss ~0.1 mm at 3e-3 rad off 180 deg,
    # ~0.1 m at 3e-4, ~10 m at 3e-5), so refuse rather than return a quietly
    # wrong answer.
    if abs(np.sin(dth)) < min_sin:
        raise ValueError("transfer angle within %.3g rad of 0/180 deg" % min_sin)

    A = np.sin(dth) * np.sqrt(n1 * n2 / (1.0 - cos_dth))

    def y(z):
        return n1 + n2 + A * (z * stumpff_s(z) - 1.0) / np.sqrt(stumpff_c(z))

    def F(z):
        yz = y(z)
        return (yz / stumpff_c(z)) ** 1.5 * stumpff_s(z) + A * np.sqrt(yz) - np.sqrt(mu) * tof

    # Upper bound: F -> +inf as z -> (2 pi)^2 (single-rev ellipse limit).
    z_hi = 4.0 * np.pi**2 - 1e-6

    # Lower bound: walk left until F < 0. For A > 0 the branch ends where
    # y = 0, and F there equals -sqrt(mu) tof < 0, so stop just inside it.
    z_lo = -4.0 * np.pi**2
    while True:
        if y(z_lo) <= 0:
            z0 = brentq(y, z_lo, z_hi, xtol=1e-14)
            z_lo = z0 + 1e-9 * max(1.0, abs(z0))
            break
        if F(z_lo) < 0:
            break
        z_lo *= 2.0
        if z_lo < -1e8:
            raise RuntimeError("no hyperbolic bracket; tof too short?")

    z = brentq(F, z_lo, z_hi, xtol=1e-14, rtol=1e-14, maxiter=500)

    yz = y(z)
    f = 1.0 - yz / n1
    g = A * np.sqrt(yz / mu)
    gdot = 1.0 - yz / n2
    v1 = (r2 - f * r1) / g
    v2 = (gdot * r2 - r1) / g
    return v1, v2


def propagate(r0, v0, t, mu=MU_EARTH):
    """Two-body numerical propagation (DOP853, tight tolerances)."""
    def rhs(_, s):
        r = s[:3]
        return np.concatenate((s[3:], -mu * r / np.linalg.norm(r) ** 3))

    sol = solve_ivp(rhs, (0.0, t), np.concatenate((r0, v0)),
                    method="DOP853", rtol=1e-12, atol=1e-12)
    return sol.y[:3, -1], sol.y[3:, -1]
