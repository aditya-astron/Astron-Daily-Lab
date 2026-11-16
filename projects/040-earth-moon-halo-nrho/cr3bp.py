# Astron Built this AP
# Project 040 | Orbital mechanics | Earth-Moon L2 halo family to the 9:2 NRHO
# Author: Aditya (aditya-astron)
# Date: 2026-11-16 IST
"""
Circular restricted three-body problem in the barycentric rotating frame.

Units: distance = Earth-Moon separation, time = 1/n (n = mean motion),
mass = M_earth + M_moon. The Earth sits at (-mu, 0, 0), the Moon at
(1 - mu, 0, 0), and the frame rotates at unit rate about +z.

Equations of motion (U = pseudo-potential):
    x'' - 2 y' = U_x,   y'' + 2 x' = U_y,   z'' = U_z
    U = (x^2 + y^2)/2 + (1 - mu)/r1 + mu/r2
Jacobi constant C = 2U - v^2 is the only integral.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# Earth-Moon system (GM values in km^3/s^2, mean distance in km)
GM_EARTH = 398600.435
GM_MOON = 4902.800
MU = GM_MOON / (GM_EARTH + GM_MOON)       # 0.0121505...
L_STAR = 384400.0                          # km
T_STAR = np.sqrt(L_STAR**3 / (GM_EARTH + GM_MOON))  # s
R_MOON = 1737.4                            # km
SYNODIC_MONTH = 29.530589 * 86400.0        # s

RTOL = 1e-12
ATOL = 1e-12


def libration_points(mu=MU):
    """Return gamma1, gamma2 (distance from the Moon) and x of L1, L2.

    From the collinear equilibrium condition U_x = 0 written as quintics
    in gamma (Szebehely):
        L1: g^5 - (3-mu) g^4 + (3-2mu) g^3 - mu g^2 + 2mu g - mu = 0
        L2: g^5 + (3-mu) g^4 + (3-2mu) g^3 - mu g^2 - 2mu g - mu = 0
    """
    q1 = lambda g: g**5 - (3 - mu) * g**4 + (3 - 2 * mu) * g**3 - mu * g**2 + 2 * mu * g - mu
    q2 = lambda g: g**5 + (3 - mu) * g**4 + (3 - 2 * mu) * g**3 - mu * g**2 - 2 * mu * g - mu
    g1 = brentq(q1, 1e-6, 0.5, xtol=1e-15)
    g2 = brentq(q2, 1e-6, 0.5, xtol=1e-15)
    return g1, g2, 1 - mu - g1, 1 - mu + g2


def potential_gradient(r, mu=MU):
    x, y, z = r
    r1 = np.sqrt((x + mu) ** 2 + y**2 + z**2)
    r2 = np.sqrt((x - 1 + mu) ** 2 + y**2 + z**2)
    k1 = (1 - mu) / r1**3
    k2 = mu / r2**3
    return np.array([x - k1 * (x + mu) - k2 * (x - 1 + mu),
                     y - k1 * y - k2 * y,
                     -k1 * z - k2 * z])


def potential_hessian(r, mu=MU):
    x, y, z = r
    d1 = np.array([x + mu, y, z])
    d2 = np.array([x - 1 + mu, y, z])
    r1 = np.linalg.norm(d1)
    r2 = np.linalg.norm(d2)
    H = np.diag([1.0, 1.0, 0.0])
    H -= (1 - mu) * (np.eye(3) / r1**3 - 3 * np.outer(d1, d1) / r1**5)
    H -= mu * (np.eye(3) / r2**3 - 3 * np.outer(d2, d2) / r2**5)
    return H


def jacobi(s, mu=MU):
    x, y, z, vx, vy, vz = s
    r1 = np.sqrt((x + mu) ** 2 + y**2 + z**2)
    r2 = np.sqrt((x - 1 + mu) ** 2 + y**2 + z**2)
    U = 0.5 * (x**2 + y**2) + (1 - mu) / r1 + mu / r2
    return 2 * U - (vx**2 + vy**2 + vz**2)


def eom(t, s, mu=MU):
    g = potential_gradient(s[:3], mu)
    return np.array([s[3], s[4], s[5],
                     2 * s[4] + g[0], -2 * s[3] + g[1], g[2]])


def jacobian_A(r, mu=MU):
    """Linearised dynamics matrix A(t) with Phi' = A Phi."""
    A = np.zeros((6, 6))
    A[:3, 3:] = np.eye(3)
    A[3:, :3] = potential_hessian(r, mu)
    A[3, 4] = 2.0
    A[4, 3] = -2.0
    return A


def eom_stm(t, w, mu=MU):
    s = w[:6]
    Phi = w[6:].reshape(6, 6)
    return np.concatenate([eom(t, s, mu), (jacobian_A(s[:3], mu) @ Phi).ravel()])


def propagate(s0, tf, mu=MU, stm=False, dense=False):
    if stm:
        w0 = np.concatenate([s0, np.eye(6).ravel()])
        sol = solve_ivp(eom_stm, (0, tf), w0, method="DOP853", rtol=RTOL,
                        atol=ATOL, args=(mu,), dense_output=dense)
        return sol.y[:6, -1], sol.y[6:, -1].reshape(6, 6), sol
    sol = solve_ivp(eom, (0, tf), s0, method="DOP853", rtol=RTOL, atol=ATOL,
                    args=(mu,), dense_output=dense)
    return sol.y[:, -1], sol


def half_period_crossing(s0, mu=MU, t_max=10.0):
    """Propagate s0 (on y = 0) with its STM to the next x-z plane crossing.

    The crossing direction is opposite to the starting vy, which skips the
    trivial root at t = 0. Returns (t_half, state, Phi).
    """
    def ev(t, w, mu):
        return w[1]
    ev.terminal = True
    ev.direction = -np.sign(s0[4])
    w0 = np.concatenate([s0, np.eye(6).ravel()])
    sol = solve_ivp(eom_stm, (0, t_max), w0, method="DOP853", rtol=RTOL,
                    atol=ATOL, args=(mu,), events=ev)
    if sol.t_events[0].size == 0:
        raise RuntimeError("no x-z plane crossing before t_max")
    w = sol.y_events[0][0]
    return sol.t_events[0][0], w[:6], w[6:].reshape(6, 6)


def symplectic_form():
    """Constant form Omega with Phi^T Omega Phi = Omega for (r, v) states.

    The canonical momenta are p = v + K r with K = [[0,-1,0],[1,0,0],[0,0,0]]
    (p_x = vx - y, p_y = vy + x), so the STM in (r, v) is symplectic with
    respect to T^T J0 T where T maps (r, v) to (q, p).
    """
    K = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 0.0]])
    I3 = np.eye(3)
    Z = np.zeros((3, 3))
    T = np.block([[I3, Z], [K, I3]])
    J0 = np.block([[Z, I3], [-I3, Z]])
    return T.T @ J0 @ T


def monodromy_from_half(Phi_half):
    """Full-period monodromy from the half-period STM.

    The orbits here are symmetric under G = diag(1,-1,1,-1,1,-1) with time
    reversal, which gives Phi(T) = G Phi(T/2)^-1 G Phi(T/2). The inverse
    uses the symplectic identity Phi^-1 = Omega^-1 Phi^T Omega, which is
    exact and avoids a numerical matrix inversion.
    """
    G = np.diag([1, -1, 1, -1, 1, -1.0])
    Om = symplectic_form()
    Phi_inv = np.linalg.solve(Om, Phi_half.T @ Om)
    return G @ Phi_inv @ G @ Phi_half
