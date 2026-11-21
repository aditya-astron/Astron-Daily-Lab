# Astron Built this AP
# Project 054 | Attitude control | Four-wheel pyramid momentum envelope and zero-crossing avoidance
# Author: Aditya (aditya-astron)
# Date: 2026-11-21 IST
"""
Momentum envelope and wheel-speed distribution for a four-wheel pyramid.

Geometry: wheel i spins about
    a_i = [cos(b) cos(phi_i), cos(b) sin(phi_i), sin(b)],  phi_i = 45, 135, 225, 315 deg
so the body momentum is h = A x, A = [a_1 .. a_4] (3x4), x_i the wheel momenta,
|x_i| <= hmax. The reachable set is the zonotope Z = A [-hmax, hmax]^4.

Envelope radius along a unit direction d is computed two independent ways:
  * LP:     max t  s.t.  A x = t d,  |x_i| <= hmax
  * facets: Z = intersection of half-spaces n.h <= hZ(n) over facet normals
            n = +-(a_i x a_j), with support hZ(n) = hmax sum_k |a_k . n|, so
            r(d) = min_{n.d > 0} hZ(n) / (n.d)
The largest sphere inside Z (guaranteed capacity in any direction) is the
distance to the nearest facet plane, min_n hZ(n)/|n|.

The 4th wheel gives a one-dimensional null space n = [1, -1, 1, -1]
(A n = 0 for any b). Wheel momenta are w = A^+ h + lam n, and lam is free:
it moves the wheels without disturbing the body. Choosing lam so that every
n_i w_i stays in [hmin, hmax] keeps all wheels away from zero speed. With
v_i = n_i w_i the reachable set is A diag(n) [hmin, hmax]^4; its centre
A n (hmin + hmax)/2 is zero, so it is exactly Z scaled by
(hmax - hmin) / (2 hmax). That is the price of never crossing zero.
"""
from itertools import combinations

import numpy as np
from scipy.optimize import linprog, minimize_scalar

PHI = np.deg2rad([45.0, 135.0, 225.0, 315.0])
NULL = np.array([1.0, -1.0, 1.0, -1.0])


def pyramid(beta):
    """3x4 matrix of wheel spin axes for elevation beta [rad]."""
    cb, sb = np.cos(beta), np.sin(beta)
    return np.vstack([cb * np.cos(PHI), cb * np.sin(PHI), sb * np.ones(4)])


def facet_normals(A):
    """Both orientations of every pairwise cross product (zonotope facet normals)."""
    ns = []
    for i, j in combinations(range(A.shape[1]), 2):
        n = np.cross(A[:, i], A[:, j])
        ns.append(n)
        ns.append(-n)
    return np.array(ns)


def support(A, n, hmax=1.0):
    return hmax * np.abs(A.T @ n).sum(axis=0)


def radius_facets(A, d, hmax=1.0):
    ns = facet_normals(A)
    nd = ns @ d
    ok = nd > 1e-12
    return np.min(support(A, ns[ok].T, hmax) / nd[ok])


def radius_lp(A, d, lo=None, hi=None):
    """Max t with A x = t d and lo <= x <= hi (elementwise bounds)."""
    m = A.shape[1]
    lo = -np.ones(m) if lo is None else lo
    hi = np.ones(m) if hi is None else hi
    c = np.zeros(m + 1)
    c[-1] = -1.0
    Aeq = np.hstack([A, -d.reshape(3, 1)])
    bounds = list(zip(lo, hi)) + [(0, None)]
    res = linprog(c, A_eq=Aeq, b_eq=np.zeros(3), bounds=bounds, method="highs")
    if res.status != 0:
        raise RuntimeError(res.message)
    return res.x[-1]


def sphere_radius(A, hmax=1.0):
    """Radius of the largest origin-centred sphere inside the envelope."""
    ns = facet_normals(A)
    return np.min(support(A, ns.T, hmax) / np.linalg.norm(ns, axis=1))


def pinv_radius(A, d, hmax=1.0):
    """Reach along d when wheels are commanded by the minimum-norm law alone."""
    return hmax / np.max(np.abs(np.linalg.pinv(A) @ d))


def best_beta():
    """Elevation that maximises the guaranteed spherical capacity."""
    r = minimize_scalar(lambda b: -sphere_radius(pyramid(b)),
                        bounds=(np.deg2rad(5), np.deg2rad(85)), method="bounded",
                        options={"xatol": 1e-10})
    return r.x, -r.fun


def fibonacci_sphere(n):
    k = np.arange(n) + 0.5
    z = 1 - 2 * k / n
    t = np.pi * (1 + 5 ** 0.5) * k
    s = np.sqrt(1 - z * z)
    return np.column_stack([s * np.cos(t), s * np.sin(t), z])


# ---------------------------------------------------------------- null motion
def bias_target(p, hmin, hmax):
    """
    Null coordinate that centres every n_i w_i inside [hmin, hmax].

    n_i (p_i + lam n_i) = n_i p_i + lam, so the bounds on lam are
    [hmin - min(n p), hmax - max(n p)]; the midpoint maximises the margin to
    both. Feasible iff range(n p) <= hmax - hmin.
    """
    q = NULL * p
    lo, hi = hmin - q.min(), hmax - q.max()
    return 0.5 * (lo + hi), lo <= hi


def simulate(A, h_of_t, t_end, dt, hmin, hmax, tau_null, k_null, use_null=True):
    """
    Wheel momenta along a body-momentum history.

    p(t) = A^+ h(t) is the minimum-norm part. lam follows the bias target
    through a saturated first-order law, lam' = sat(k (lam* - lam), tau_null),
    which is the null torque each wheel spends on redistribution.
    """
    Ap = np.linalg.pinv(A)
    t = np.arange(0.0, t_end + dt / 2, dt)
    W = np.empty((t.size, 4))
    H = np.empty((t.size, 3))
    lam = 0.0
    feasible = True
    for k, tk in enumerate(t):
        h = h_of_t(tk)
        p = Ap @ h
        W[k] = p + lam * NULL
        H[k] = h
        if use_null:
            target, ok = bias_target(p, hmin, hmax)
            feasible &= ok
            lam += dt * np.clip(k_null * (target - lam), -tau_null, tau_null)
    return t, W, H, feasible


def zero_crossings(W):
    s = np.sign(W)
    return int(np.sum(s[1:] * s[:-1] < 0) + np.sum(s[1:] == 0))
