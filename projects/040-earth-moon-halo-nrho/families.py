# Astron Built this AP
# Project 040 | Orbital mechanics | Earth-Moon L2 halo family to the 9:2 NRHO
# Author: Aditya (aditya-astron)
# Date: 2026-11-16 IST
"""
Differential correction and continuation of symmetric periodic orbits.

Every orbit here is symmetric about the x-z plane, so it is fully defined
by a perpendicular crossing  s0 = (x0, 0, z0, 0, vy0, 0)  and is periodic
when the next crossing is also perpendicular: vx(T/2) = vz(T/2) = 0.

Single-shooting Jacobian. The crossing time floats with the initial state,
so the partials of the stopped state include the time shift dt = -dy/vy:
    d s(T/2) / d X = Phi[:, X] - f(s) / vy * Phi[1, X]
where f is the vector field and X picks the free initial components.
"""

from dataclasses import dataclass

import numpy as np

from cr3bp import MU, eom, half_period_crossing, monodromy_from_half, jacobi

IX, IZ, IVY = 0, 2, 4


@dataclass
class Orbit:
    s0: np.ndarray       # perpendicular x-z crossing state
    period: float
    Phi_half: np.ndarray
    s_half: np.ndarray

    @property
    def C(self):
        return jacobi(self.s0)

    def monodromy(self):
        return monodromy_from_half(self.Phi_half)


def _crossing_partials(s0, mu):
    th, sh, Ph = half_period_crossing(s0, mu)
    f = eom(0.0, sh, mu)
    D = Ph - np.outer(f, Ph[1, :]) / sh[4]
    return th, sh, Ph, D


def _assemble(X, s_template, free):
    s = s_template.copy()
    s[free] = X
    return s


def correct(s0, free=(IX, IVY), targets=(3, 5), extra=None, mu=MU,
            tol=1e-11, max_iter=30):
    """Newton solve with the given free initial components.

    targets: components of the half-period state driven to zero (vx, vz).
    extra: optional (residual(X, th), gradient(X, th, dth)) pair added as one
           more equation, used for pseudo-arclength and period fixing; dth is
           the gradient of the half-period with respect to X.
    The system is square when len(free) == len(targets) + (extra is not None).
    """
    free = list(free)
    targets = list(targets)
    X = np.asarray(s0, float)[free].copy()
    s_t = np.asarray(s0, float).copy()
    for _ in range(max_iter):
        s = _assemble(X, s_t, free)
        th, sh, Ph, D = _crossing_partials(s, mu)
        F = [sh[i] for i in targets]
        J = [D[i, free] for i in targets]
        # The half-period crossing time also moves with X.
        dth = -Ph[1, free] / sh[4]
        if extra is not None:
            res, grad = extra
            F.append(res(X, th))
            J.append(grad(X, th, dth))
        F = np.array(F)
        J = np.array(J)
        if np.max(np.abs(F)) < tol:
            return Orbit(s, 2 * th, Ph, sh)
        X = X - np.linalg.lstsq(J, F, rcond=None)[0]
    raise RuntimeError(f"differential correction did not converge, |F| = {np.max(np.abs(F)):.2e}")


def correct_fixed_period(s0, T, mu=MU, tol=1e-11):
    """Halo orbit with prescribed period, free (x0, z0, vy0)."""
    extra = (lambda X, th: th - T / 2,
             lambda X, th, dth: dth)
    return correct(s0, free=(IX, IZ, IVY), extra=extra, mu=mu, tol=tol)


def family_tangent(orb, free, targets=(3, 5), mu=MU):
    """Null vector of the constraint Jacobian at a converged orbit."""
    _, _, _, D = _crossing_partials(orb.s0, mu)
    J = D[np.ix_(list(targets), list(free))]
    _, _, Vt = np.linalg.svd(J)
    return Vt[-1]


def continue_family(orb, ds, n_steps, free=(IX, IZ, IVY), targets=(3, 5),
                    direction=None, stop=None, mu=MU):
    """Pseudo-arclength continuation.

    Each step predicts X* + ds n and corrects with the extra equation
    (X - X*) . n - ds = 0, so the solver can follow folds where a natural
    parameter (such as Az) turns back.
    """
    free = list(free)
    fam = [orb]
    n_prev = family_tangent(orb, free, targets, mu)
    if direction is not None and np.dot(n_prev, direction) < 0:
        n_prev = -n_prev
    for _ in range(n_steps):
        cur = fam[-1]
        n = family_tangent(cur, free, targets, mu)
        if np.dot(n, n_prev) < 0:
            n = -n
        X_star = cur.s0[free].copy()
        guess = cur.s0.copy()
        guess[free] = X_star + ds * n
        extra = (lambda X, th, Xs=X_star, nn=n: np.dot(X - Xs, nn) - ds,
                 lambda X, th, dth, nn=n: nn)
        nxt = correct(guess, free=free, targets=targets, extra=extra, mu=mu)
        fam.append(nxt)
        n_prev = n
        if stop is not None and stop(nxt):
            break
    return fam


def stability(orb):
    """Monodromy eigenvalues and the two Broucke stability indices.

    Eigenvalues come in reciprocal pairs and include a trivial pair at 1.
    Each nontrivial pair gives nu = (lam + 1/lam)/2; |nu| <= 1 is
    (linearly) stable in that mode.
    """
    M = orb.monodromy()
    ev = np.linalg.eigvals(M)
    # Remove the trivial pair: the two eigenvalues closest to 1.
    idx = np.argsort(np.abs(ev - 1.0))
    rest = ev[idx[2:]]
    nus = []
    used = np.zeros(4, bool)
    for i in range(4):
        if used[i]:
            continue
        j = min((j for j in range(4) if j != i and not used[j]),
                key=lambda j: abs(rest[i] * rest[j] - 1))
        used[i] = used[j] = True
        nus.append(((rest[i] + rest[j]) / 2).real)
    return ev, sorted(nus, key=abs, reverse=True)


def lyapunov_out_of_plane_trace(orb):
    """Trace of the (z, vz) block of a planar orbit's monodromy.

    For a planar orbit the z-motion decouples, so this 2x2 block has det 1
    and its eigenvalues reach +1 exactly when the trace equals 2. That is
    where the halo family branches off the Lyapunov family.
    """
    M = orb.monodromy()
    return M[2, 2] + M[5, 5]
