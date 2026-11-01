# Astron Built this AP
# Project 004 | Spacecraft power systems | Solar array MPPT through a LEO eclipse cycle
# Author: Aditya (aditya-astron)
# Date: 2026-11-01 IST
"""Single-diode model of a triple-junction GaAs cell, with parameters extracted
from the four numbers a datasheet gives (Isc, Voc, Vmp, Imp), plus the array
built from it.

    I = Iph - I0 (exp((V + I Rs)/a) - 1) - (V + I Rs)/Rsh,   a = n k T / q

The equation is implicit in I. It has an exact explicit form through the
Lambert W function; the argument of W overflows a double for any realistic
cell, so it is evaluated as the Wright omega function, omega(x) = W(exp(x)).
"""

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq
from scipy.special import wrightomega

K_B = 1.380649e-23    # J/K
Q_E = 1.602176634e-19  # C
T0_K = 273.15
G_AM0 = 1367.0        # W/m^2, reference irradiance for the cell data below


@dataclass(frozen=True)
class CellDatasheet:
    """Beginning-of-life AM0 numbers at 28 C. Representative of a commercial
    ~29.5 % triple-junction GaAs cell (30 cm^2 class), not a specific part."""
    isc: float = 0.5196     # A
    voc: float = 2.700      # V
    vmp: float = 2.411      # V
    imp: float = 0.5029     # A
    t_ref_c: float = 28.0
    area_m2: float = 30.18e-4
    di_sc_dt: float = 0.36e-3   # A/K
    dv_oc_dt: float = -6.2e-3   # V/K
    # Effective ideality of three junctions in series. With this fill factor
    # (0.864) the fit has Rsh > 0 only for n below ~3.05; 2.9 gives
    # Rs ~ 53 mOhm and Rsh ~ 2.8 kOhm, both typical for a cell this size.
    n_ideality: float = 2.9


@dataclass(frozen=True)
class CellParams:
    iph: float
    i0: float
    rs: float
    rsh: float
    a: float   # modified thermal voltage n k T / q, V


def thermal_voltage(n, t_c):
    return n * K_B * (t_c + T0_K) / Q_E


def cell_current(v, p: CellParams):
    """Explicit single-diode solution, vectorised over v."""
    v = np.asarray(v, dtype=float)
    rs, rsh, a = p.rs, p.rsh, p.a
    rsum = rs + rsh
    x = (math.log(rs * p.i0 * rsh / (a * rsum))
         + rsh * (rs * (p.iph + p.i0) + v) / (a * rsum))
    w = np.real(wrightomega(x))
    return (rsh * (p.iph + p.i0) - v) / rsum - (a / rs) * w


def cell_current_newton(v, p: CellParams, iters=60):
    """Newton on the implicit equation. Used only to cross-check the W form."""
    i = p.iph
    for _ in range(iters):
        vd = v + i * p.rs
        e = math.exp(vd / p.a)
        f = p.iph - p.i0 * (e - 1.0) - vd / p.rsh - i
        df = -p.i0 * e * p.rs / p.a - p.rs / p.rsh - 1.0
        step = f / df
        i -= step
        if abs(step) < 1e-15:
            break
    return i


def extract_params(ds: CellDatasheet = CellDatasheet()) -> CellParams:
    """Solve for Iph, I0, Rs, Rsh with n fixed, from four conditions:
    I(0) = Isc, I(Voc) = 0, I(Vmp) = Imp, and dP/dV = 0 at (Vmp, Imp).

    For a trial Rs the first three are linear in (Iph, I0, 1/Rsh):
        Iph - I0 (exp(Vd/a) - 1) - Vd / Rsh = I,   Vd = V + I Rs
    so they are solved exactly as a 3x3 system. What remains is the
    maximum-power condition. Differentiating the implicit equation gives
        dI/dV = -G / (1 + G Rs),   G = (I0/a) exp(Vd/a) + 1/Rsh
    and dP/dV = 0 means G / (1 + G Rs) = Imp / Vmp. That is one equation in
    Rs, solved with a bracketed Brent root.
    """
    a = thermal_voltage(ds.n_ideality, ds.t_ref_c)
    pts = [(0.0, ds.isc), (ds.voc, 0.0), (ds.vmp, ds.imp)]

    def linear_part(rs):
        m = np.empty((3, 3))
        rhs = np.empty(3)
        for k, (v, i) in enumerate(pts):
            vd = v + i * rs
            m[k] = [1.0, -math.expm1(vd / a), -vd]
            rhs[k] = i
        iph, i0, gsh = np.linalg.solve(m, rhs)
        return iph, i0, gsh

    def mpp_residual(rs):
        iph, i0, gsh = linear_part(rs)
        g = i0 / a * math.exp((ds.vmp + ds.imp * rs) / a) + gsh
        return g / (1.0 + g * rs) * ds.vmp / ds.imp - 1.0

    # Scan for the bracket where the fit still has I0 > 0 and Rsh > 0.
    grid = np.linspace(1e-6, (ds.voc - ds.vmp) / ds.imp, 400)
    for lo, hi in zip(grid[:-1], grid[1:]):
        try:
            change = mpp_residual(lo) * mpp_residual(hi) < 0
        except np.linalg.LinAlgError:
            continue
        if change:
            rs = brentq(mpp_residual, lo, hi, xtol=1e-15)
            break
    else:
        raise RuntimeError("no Rs satisfies dP/dV = 0 at the datasheet MPP")
    iph, i0, gsh = linear_part(rs)
    if i0 <= 0 or gsh <= 0:
        raise RuntimeError(f"non-physical fit: I0={i0}, 1/Rsh={gsh} for n={ds.n_ideality}")
    return CellParams(iph, i0, rs, 1.0 / gsh, a)


def params_at(ref: CellParams, ds: CellDatasheet, g_wm2, t_c) -> CellParams:
    """Move the reference fit to another irradiance and temperature.

    Iph scales with irradiance and follows the Isc coefficient. I0 is then
    re-derived so that Voc at full sun follows the datasheet dVoc/dT exactly
    (the Villalva approach), which sidesteps the band-gap physics of three
    stacked junctions. Rs and Rsh are held constant.
    """
    dt = t_c - ds.t_ref_c
    a = thermal_voltage(ds.n_ideality, t_c)
    # Iph at reference is slightly above Isc because of the Rs/Rsh drop.
    iph_full = ref.iph * (1.0 + ds.di_sc_dt * dt / ds.isc)
    voc_full = ds.voc + ds.dv_oc_dt * dt
    i0 = (iph_full - voc_full / ref.rsh) / (math.exp(voc_full / a) - 1.0)
    return CellParams(iph_full * g_wm2 / G_AM0, i0, ref.rs, ref.rsh, a)


@dataclass(frozen=True)
class Array:
    """n_series cells per string, n_parallel identical strings, no mismatch.
    Blocking-diode and harness drops are ignored."""
    n_series: int = 16
    n_parallel: int = 12

    def current(self, v, p: CellParams):
        return self.n_parallel * cell_current(np.asarray(v) / self.n_series, p)

    def power(self, v, p: CellParams):
        v = np.asarray(v, dtype=float)
        return v * self.current(v, p)

    def voc(self, p: CellParams):
        return self.n_series * open_circuit_voltage(p)


def open_circuit_voltage(p: CellParams):
    """Voc from I = 0: Iph = I0 (exp(V/a) - 1) + V/Rsh, solved by Newton."""
    v = p.a * math.log(p.iph / p.i0 + 1.0)
    for _ in range(50):
        e = math.exp(v / p.a)
        f = p.i0 * (e - 1.0) + v / p.rsh - p.iph
        df = p.i0 * e / p.a + 1.0 / p.rsh
        step = f / df
        v -= step
        if abs(step) < 1e-14:
            break
    return v


def true_mpp(array: Array, p_list, iters=80):
    """Vectorised golden-section search for the maximum power point of each
    operating condition in p_list. P(V) is unimodal on [0, Voc] for a
    uniformly lit single-diode array, so golden section is enough."""
    n = len(p_list)
    iph = np.array([p.iph for p in p_list])
    i0 = np.array([p.i0 for p in p_list])
    a = np.array([p.a for p in p_list])
    rs = p_list[0].rs
    rsh = p_list[0].rsh
    rsum = rs + rsh
    ln_k = np.log(rs * i0 * rsh / (a * rsum))

    def power(v):
        vc = v / array.n_series
        x = ln_k + rsh * (rs * (iph + i0) + vc) / (a * rsum)
        ic = (rsh * (iph + i0) - vc) / rsum - (a / rs) * np.real(wrightomega(x))
        return v * array.n_parallel * ic

    voc = np.array([array.voc(p) if p.iph > 0 else 0.0 for p in p_list])
    lo = np.zeros(n)
    hi = voc.copy()
    phi = (math.sqrt(5.0) - 1.0) / 2.0
    for _ in range(iters):
        c = hi - phi * (hi - lo)
        d = lo + phi * (hi - lo)
        left = power(c) > power(d)   # maximum lies in [lo, d]
        hi = np.where(left, d, hi)
        lo = np.where(left, lo, c)
    v = 0.5 * (lo + hi)
    pmax = np.where(voc > 0, power(v), 0.0)
    return np.where(voc > 0, v, 0.0), pmax
