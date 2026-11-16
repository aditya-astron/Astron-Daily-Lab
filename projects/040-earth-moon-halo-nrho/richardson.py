# Astron Built this AP
# Project 040 | Orbital mechanics | Earth-Moon L2 halo family to the 9:2 NRHO
# Author: Aditya (aditya-astron)
# Date: 2026-11-16 IST
"""
Richardson (1980) third-order Lindstedt-Poincare halo approximation.

Coordinates are centred on the libration point, scaled by gamma (its
distance to the Moon) and share the axes of the rotating frame. The
Legendre coefficients of the expanded potential are
    L1: c_n = [mu + (-1)^n (1-mu) g^(n+1)/(1-g)^(n+1)] / g^3
    L2: c_n = [(-1)^n mu + (-1)^n (1-mu) g^(n+1)/(1+g)^(n+1)] / g^3

The halo condition comes from removing secular terms at third order:
    l1 Ax^2 + l2 Az^2 + Delta = 0,   Delta = lambda^2 - c2
which sets the minimum in-plane amplitude sqrt(-Delta/l1) at Az = 0.
"""

from dataclasses import dataclass

import numpy as np

from cr3bp import MU, libration_points


@dataclass
class RichardsonCoeffs:
    gamma: float
    x_lib: float
    c2: float
    c3: float
    c4: float
    lam: float
    k: float
    delta: float
    a21: float
    a22: float
    a23: float
    a24: float
    a31: float
    a32: float
    b21: float
    b22: float
    b31: float
    b32: float
    d21: float
    d31: float
    d32: float
    s1: float
    s2: float
    l1: float
    l2: float


def coefficients(point=2, mu=MU, gamma=None, x_lib=None):
    if gamma is None:
        g1, g2, x1, x2 = libration_points(mu)
        gamma, x_lib = (g1, x1) if point == 1 else (g2, x2)
    g = gamma

    def c(n):
        if point == 1:
            return (mu + (-1) ** n * (1 - mu) * g ** (n + 1) / (1 - g) ** (n + 1)) / g**3
        return ((-1) ** n * mu + (-1) ** n * (1 - mu) * g ** (n + 1) / (1 + g) ** (n + 1)) / g**3

    c2, c3, c4 = c(2), c(3), c(4)
    # In-plane linear frequency: lambda^4 + (c2-2) lambda^2 - (c2-1)(1+2c2) = 0
    lam2 = (-(c2 - 2) + np.sqrt((c2 - 2) ** 2 + 4 * (c2 - 1) * (1 + 2 * c2))) / 2
    lam = np.sqrt(lam2)
    k = (lam2 + 1 + 2 * c2) / (2 * lam)
    delta = lam2 - c2

    d1 = 3 * lam2 / k * (k * (6 * lam2 - 1) - 2 * lam)
    d2 = 8 * lam2 / k * (k * (11 * lam2 - 1) - 2 * lam)

    a21 = 3 * c3 * (k**2 - 2) / (4 * (1 + 2 * c2))
    a22 = 3 * c3 / (4 * (1 + 2 * c2))
    a23 = -3 * c3 * lam / (4 * k * d1) * (3 * k**3 * lam - 6 * k * (k - lam) + 4)
    a24 = -3 * c3 * lam / (4 * k * d1) * (2 + 3 * k * lam)
    b21 = -3 * c3 * lam / (2 * d1) * (3 * k * lam - 4)
    b22 = 3 * c3 * lam / d1
    d21 = -c3 / (2 * lam2)

    a31 = (-9 * lam / (4 * d2) * (4 * c3 * (k * a23 - b21) + k * c4 * (4 + k**2))
           + (9 * lam2 + 1 - c2) / (2 * d2) * (3 * c3 * (2 * a21 - k * b21) + c4 * (2 + 3 * k**2)))
    a32 = -1 / d2 * (9 * lam / 4 * (4 * c3 * (k * a24 - b22) + k * c4)
                     + 1.5 * (9 * lam2 + 1 - c2) * (c3 * (k * b22 + d21 - 2 * a22) - c4))
    b31 = 3 / (8 * d2) * (8 * lam * (3 * c3 * (k * b21 - 2 * a23) - c4 * (2 + 3 * k**2))
                          + (9 * lam2 + 1 + 2 * c2) * (4 * c3 * (k * a23 - b21) + k * c4 * (4 + k**2)))
    b32 = 1 / d2 * (9 * lam * (c3 * (k * b22 + d21 - 2 * a22) - c4)
                    + 3 / 8 * (9 * lam2 + 1 + 2 * c2) * (4 * c3 * (k * a24 - b22) + k * c4))
    d31 = 3 / (64 * lam2) * (4 * c3 * a24 + c4)
    d32 = 3 / (64 * lam2) * (4 * c3 * (a23 - d21) + c4 * (4 + k**2))

    den = 2 * lam * (lam * (1 + k**2) - 2 * k)
    s1 = (1.5 * c3 * (2 * a21 * (k**2 - 2) - a23 * (k**2 + 2) - 2 * k * b21)
          - 3 / 8 * c4 * (3 * k**4 - 8 * k**2 + 8)) / den
    s2 = (1.5 * c3 * (2 * a22 * (k**2 - 2) + a24 * (k**2 + 2) + 2 * k * b22 + 5 * d21)
          + 3 / 8 * c4 * (12 - k**2)) / den
    a1 = -1.5 * c3 * (2 * a21 + a23 + 5 * d21) - 3 / 8 * c4 * (12 - k**2)
    a2 = 1.5 * c3 * (a24 - 2 * a22) + 9 / 8 * c4
    l1 = a1 + 2 * lam2 * s1
    l2 = a2 + 2 * lam2 * s2

    return RichardsonCoeffs(g, x_lib, c2, c3, c4, lam, k, delta, a21, a22, a23,
                            a24, a31, a32, b21, b22, b31, b32, d21, d31, d32,
                            s1, s2, l1, l2)


def ax_from_az(rc, Az):
    """In-plane amplitude (gamma units) that makes the orbit periodic."""
    val = (-rc.l2 * Az**2 - rc.delta) / rc.l1
    if val < 0:
        raise ValueError("no halo orbit for this Az in third-order theory")
    return np.sqrt(val)


def halo_state(rc, Az_km, L_star, tau1, south=True):
    """Rotating-frame state on the third-order halo with out-of-plane
    amplitude Az_km at phase tau1. Returns (state, period).

    'south' means the far (apolune) side lies below the Moon's orbit plane,
    the branch that bends into the Gateway-type L2 NRHO.
    """
    Az = Az_km / (rc.gamma * L_star)
    Ax = ax_from_az(rc, Az)
    # delta_m = +1 puts z = +Az at tau1 = 0 (the Moon side), so the far
    # side at tau1 = pi has z < 0.
    dm = 1.0 if south else -1.0
    omega = 1 + rc.s1 * Ax**2 + rc.s2 * Az**2
    c1, c2_, c3_ = np.cos(tau1), np.cos(2 * tau1), np.cos(3 * tau1)
    s1_, s2_, s3_ = np.sin(tau1), np.sin(2 * tau1), np.sin(3 * tau1)

    x = (rc.a21 * Ax**2 + rc.a22 * Az**2 - Ax * c1
         + (rc.a23 * Ax**2 - rc.a24 * Az**2) * c2_
         + (rc.a31 * Ax**3 - rc.a32 * Ax * Az**2) * c3_)
    y = (rc.k * Ax * s1_ + (rc.b21 * Ax**2 - rc.b22 * Az**2) * s2_
         + (rc.b31 * Ax**3 - rc.b32 * Ax * Az**2) * s3_)
    z = (dm * Az * c1 + dm * rc.d21 * Ax * Az * (c2_ - 3)
         + dm * (rc.d32 * Az * Ax**2 - rc.d31 * Az**3) * c3_)

    # d/dt = lambda * omega * d/dtau1 (time in rotating-frame units)
    w = rc.lam * omega
    vx = w * (Ax * s1_ - 2 * (rc.a23 * Ax**2 - rc.a24 * Az**2) * s2_
              - 3 * (rc.a31 * Ax**3 - rc.a32 * Ax * Az**2) * s3_)
    vy = w * (rc.k * Ax * c1 + 2 * (rc.b21 * Ax**2 - rc.b22 * Az**2) * c2_
              + 3 * (rc.b31 * Ax**3 - rc.b32 * Ax * Az**2) * c3_)
    vz = w * (-dm * Az * s1_ - 2 * dm * rc.d21 * Ax * Az * s2_
              - 3 * dm * (rc.d32 * Az * Ax**2 - rc.d31 * Az**3) * s3_)

    g = rc.gamma
    state = np.array([rc.x_lib + g * x, g * y, g * z, g * vx, g * vy, g * vz])
    period = 2 * np.pi / w
    return state, period
