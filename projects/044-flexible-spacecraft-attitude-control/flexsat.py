# Astron Built this AP
# Project 044 | PID and control systems | Attitude control of a spacecraft with flexible solar arrays
# Author: Aditya (aditya-astron)
# Date: 2026-11-17 IST
"""
Single-axis structural model: a rigid hub of radius r carrying two identical
cantilevered panels (Euler-Bernoulli beams, length l, mass/length m, stiffness
EI) that bend in the plane of rotation. Only the antisymmetric bending couples
to hub rotation, so one panel's deflection is mirrored onto the other.

Displacement of a panel point at distance x from its root:
    y(x, t) = (r + x) theta(t) + w(x, t),    w(0) = w'(0) = 0
w is the elastic deflection in the hub frame (clamped at the root).

Hybrid equations (hub-fixed / constrained modes w = sum phi_k eta_k):
    J_t theta'' + sum_k L_k eta_k''                      = tau
    L_k theta'' + eta_k'' + 2 z w_k eta_k' + w_k^2 eta_k = 0
J_t is the whole-body inertia and L_k = 2 int m (r + x) phi_k dx the modal
rotational participation of mode k (both panels). With one mode the hub
transfer function theta/tau has zeros at the constrained (cantilever)
frequency w_c and poles at the free-free frequency
    w_p = w_c / sqrt(1 - L^2 / J_t)
"""

import numpy as np
from scipy.linalg import eigh

BETA1 = 1.8751040687119611   # first root of cos(b) cosh(b) = -1


def beam_element(le, m, EI):
    """Hermite cubic beam element mass and stiffness, DOFs (w1, w1', w2, w2')."""
    k = EI / le**3 * np.array([[12, 6 * le, -12, 6 * le],
                               [6 * le, 4 * le**2, -6 * le, 2 * le**2],
                               [-12, -6 * le, 12, -6 * le],
                               [6 * le, 2 * le**2, -6 * le, 4 * le**2]])
    mm = m * le / 420 * np.array([[156, 22 * le, 54, -13 * le],
                                  [22 * le, 4 * le**2, 13 * le, -3 * le**2],
                                  [54, 13 * le, 156, -22 * le],
                                  [-13 * le, -3 * le**2, -22 * le, 4 * le**2]])
    return mm, k


def coupling_element(le, m, x0, r):
    """int m (r + x) N_i dx over one element starting at x0 (exact: cubic*linear)."""
    g, wq = np.polynomial.legendre.leggauss(4)
    out = np.zeros(4)
    for gi, wi in zip(g, wq):
        s = 0.5 * (gi + 1.0)
        n = np.array([1 - 3 * s**2 + 2 * s**3, le * (s - 2 * s**2 + s**3),
                      3 * s**2 - 2 * s**3, le * (-s**2 + s**3)])
        out += 0.5 * wi * m * (r + x0 + s * le) * n * le
    return out


class FlexSat:
    """Hub + two panels, FE discretised, reduced to constrained modes."""

    def __init__(self, J_hub=600.0, r=1.0, l=6.0, m=5.0, EI=1300.0,
                 n_el=40, zeta=0.005):
        self.J_hub, self.r, self.l, self.m, self.EI = J_hub, r, l, m, EI
        self.zeta = zeta
        le = l / n_el
        nd = 2 * (n_el + 1)
        M = np.zeros((nd, nd))
        K = np.zeros((nd, nd))
        b = np.zeros(nd)
        for e in range(n_el):
            i = slice(2 * e, 2 * e + 4)
            me, ke = beam_element(le, m, EI)
            M[i, i] += me
            K[i, i] += ke
            b[i] += coupling_element(le, m, e * le, r)
        # clamp the root (w = w' = 0)
        self.Mww, self.Kww, self.b = M[2:, 2:], K[2:, 2:], b[2:]
        self.J_panel = m * ((r + l)**3 - r**3) / 3.0
        self.J_t = J_hub + 2 * self.J_panel

        # constrained modes of one panel, mass normalised
        w2, phi = eigh(self.Kww, self.Mww)
        self.w_c = np.sqrt(w2)
        self.phi = phi
        # two panels deflecting antisymmetrically: both rotate with the hub, so
        # the generalised coordinate eta_k covers both panels with mass 1:
        # phi_total = phi / sqrt(2), L_k = 2 * b.phi / sqrt(2)
        self.L = np.sqrt(2.0) * (self.b @ phi)
        self.L *= np.sign(self.L)

    def free_free_fe(self, n=3):
        """First n nonzero elastic frequencies of the full hub+panels FE model."""
        nd = self.Mww.shape[0]
        M = np.zeros((nd + 1, nd + 1))
        K = np.zeros((nd + 1, nd + 1))
        # per-panel antisymmetric half model with half the hub
        M[0, 0] = self.J_t / 2.0
        M[0, 1:] = M[1:, 0] = self.b
        M[1:, 1:] = self.Mww
        K[1:, 1:] = self.Kww
        w2 = eigh(K, M, eigvals_only=True)
        return np.sqrt(np.clip(w2[1:n + 1], 0, None))

    def modal_free_free(self, n_modes):
        """Elastic frequencies of the hybrid model truncated to n_modes."""
        L = self.L[:n_modes]
        M = np.block([[np.array([[self.J_t]]), L[None, :]],
                      [L[:, None], np.eye(n_modes)]])
        K = np.diag(np.r_[0.0, self.w_c[:n_modes]**2])
        w2 = eigh(K, M, eigvals_only=True)
        return np.sqrt(np.clip(w2[1:], 0, None))

    def state_space(self, n_modes=2, freq_scale=1.0):
        """(A, B, C) with x = [theta, eta..., theta', eta'...], outputs [theta, theta'].

        freq_scale multiplies every modal frequency (on-orbit uncertainty).
        """
        n = n_modes
        L = self.L[:n]
        wk = self.w_c[:n] * freq_scale
        M = np.block([[np.array([[self.J_t]]), L[None, :]],
                      [L[:, None], np.eye(n)]])
        K = np.diag(np.r_[0.0, wk**2])
        D = np.diag(np.r_[0.0, 2 * self.zeta * wk])
        Mi = np.linalg.inv(M)
        A = np.block([[np.zeros((n + 1, n + 1)), np.eye(n + 1)],
                      [-Mi @ K, -Mi @ D]])
        B = np.r_[np.zeros(n + 1), Mi[:, 0]][:, None]
        C = np.zeros((2, 2 * n + 2))
        C[0, 0] = 1.0
        C[1, n + 1] = 1.0
        return A, B, C

    def tip_deflection(self, x, n_modes):
        """Elastic tip deflection w(l) of one panel from state x."""
        tip = self.phi[-2, :n_modes] / np.sqrt(2.0) * np.sign(self.b @ self.phi[:, :n_modes])
        return x[..., 1:1 + n_modes] @ tip


def cantilever_mode1(x, l):
    """Analytic first clamped-free mode shape (unnormalised)."""
    b = BETA1 / l
    s = (np.cosh(BETA1) + np.cos(BETA1)) / (np.sinh(BETA1) + np.sin(BETA1))
    return np.cosh(b * x) - np.cos(b * x) - s * (np.sinh(b * x) - np.sin(b * x))
