# Astron Built this AP
# Project 050 | Rocket propulsion | Regenerative cooling of a LOX/methane thrust chamber
# Author: Aditya (aditya-astron)
# Date: 2026-11-20 IST
"""
2D finite-volume conduction in one channel/land symmetry cell.

Used to check the 1D wall + fin model at a single station. The cell spans
from the channel centreline (x = 0) to the land centreline (x = p/2), both
symmetry planes, and from the hot face (y = 0) to the closeout (y = tw + h):

    y in [0, tw]          solid across the whole width (hot-wall ligament)
    y in [tw, tw + h]     solid only for x in [w/2, p/2] (half land)

Boundary conditions
    hot face        -k dT/dy = hg (Taw - T)
    channel faces   -k dT/dn = hc (T - Tb)   (channel floor and land side)
    closeout, symmetry planes   adiabatic

Cells are uniform within each of the four sub-blocks (half channel width,
half land width, ligament, land height) so every interface lies on a face.
"""
from __future__ import annotations

import math

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def _edges(a: float, b: float, n: int) -> np.ndarray:
    return np.linspace(a, b, n + 1)


def solve_cell(w: float, b: float, tw: float, h: float, k: float,
               hg: float, Taw: float, hc: float, Tb: float, n: int = 40):
    """Return (Q per cell per unit length [W/m], hot-face T profile, x centres, T field)."""
    half_w, half_b = 0.5 * w, 0.5 * b
    if half_w <= 0.0:
        raise ValueError("channel width must be positive")
    xe = _edges(0.0, half_w, n)
    if half_b > 0.0:
        xe = np.concatenate([xe[:-1], _edges(half_w, half_w + half_b, n)])
    ye = _edges(0.0, tw, n)
    if h > 0.0:
        ye = np.concatenate([ye[:-1], _edges(tw, tw + h, n)])
    xc, yc = 0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1])
    dx, dy = np.diff(xe), np.diff(ye)
    nx, ny = len(xc), len(yc)
    solid = np.ones((nx, ny), bool)
    solid[np.ix_(xc < half_w, yc > tw)] = False

    idx = -np.ones((nx, ny), int)
    idx[solid] = np.arange(solid.sum())
    N = int(solid.sum())
    rows, cols, vals = [], [], []
    rhs = np.zeros(N)
    diag = np.zeros(N)

    def link(p, q, G):
        rows.extend([p, q]); cols.extend([q, p]); vals.extend([-G, -G])
        diag[p] += G; diag[q] += G

    def to_fluid(p, G, T):
        diag[p] += G; rhs[p] += G * T

    for i in range(nx):
        for j in range(ny):
            p = idx[i, j]
            if p < 0:
                continue
            # x-neighbour
            if i + 1 < nx:
                if solid[i + 1, j]:
                    G = k * dy[j] / (0.5 * (dx[i] + dx[i + 1]))
                    link(p, idx[i + 1, j], G)
            # left face of a land cell facing the channel
            if i > 0 and not solid[i - 1, j]:
                G = dy[j] / (0.5 * dx[i] / k + 1.0 / hc)
                to_fluid(p, G, Tb)
            # y-neighbour
            if j + 1 < ny:
                if solid[i, j + 1]:
                    G = k * dx[i] / (0.5 * (dy[j] + dy[j + 1]))
                    link(p, idx[i, j + 1], G)
                else:  # channel floor
                    G = dx[i] / (0.5 * dy[j] / k + 1.0 / hc)
                    to_fluid(p, G, Tb)
            elif h == 0.0 or xc[i] < half_w:
                # no land: top of the ligament is the channel floor everywhere
                G = dx[i] / (0.5 * dy[j] / k + 1.0 / hc)
                to_fluid(p, G, Tb)
            if j == 0:
                G = dx[i] / (0.5 * dy[0] / k + 1.0 / hg)
                to_fluid(p, G, Taw)

    A = sp.coo_matrix((vals + list(diag), (rows + list(range(N)), cols + list(range(N)))),
                      shape=(N, N)).tocsr()
    T = spla.spsolve(A, rhs)
    field = np.full((nx, ny), np.nan)
    field[solid] = T
    # heat in through the hot face
    T0 = field[:, 0]
    Gh = dx / (0.5 * dy[0] / k + 1.0 / hg)
    Q = float(np.sum(Gh * (Taw - T0)))
    # surface temperature from the half-cell resistance
    Ts = T0 + (Gh / dx) * (Taw - T0) * 0.5 * dy[0] / k
    return Q, Ts, xc, field


def fin_model_flat(w: float, b: float, tw: float, h: float, k: float,
                   hg: float, Taw: float, hc: float, Tb: float):
    """The 1D ligament + adiabatic-tip fin model on the same flat cell.

    Series resistances per unit length for a full pitch p = w + b:
        1/(hg p) + tw/(k p) + 1/(hc (w + 2 eta h)),  eta = tanh(mh)/(mh)
    Returns (Q per half cell, hot-face temperature).
    """
    p = w + b
    m = math.sqrt(2.0 * hc / (k * b))
    eta = math.tanh(m * h) / (m * h) if h > 0 else 1.0
    R = 1.0 / (hg * p) + tw / (k * p) + 1.0 / (hc * (w + 2.0 * eta * h))
    Q = (Taw - Tb) / R
    return 0.5 * Q, Taw - Q / (hg * p)
