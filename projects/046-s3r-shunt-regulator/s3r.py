# Astron Built this AP
# Project 046 | Spacecraft power systems | S3R shunt regulator for a 50 V GEO bus
# Author: Aditya (aditya-astron)
# Date: 2026-11-18 IST
"""
Sequential switching shunt regulator (S3R) building blocks.

Array side
    Each cell is a single-diode model with no shunt leakage:
        I = Iph - I0 (exp((V + I Rs)/a) - 1),   a = n k T / q
    Iph and I0 are pinned to the cell's Isc and Voc in closed form, and the
    implicit equation is solved explicitly with the Lambert W function:
        I = Iph + I0 - (a/Rs) W( (I0 Rs / a) exp((V + (Iph + I0) Rs)/a) )
    A string's harness resistance is folded into Rs per cell, and the bus
    sees the string through a blocking diode, so V_cell = (V_bus + V_d)/Ns.

Regulator side
    Section k is shorted by its shunt switch when the bus reaches hi_k and
    released when it falls to lo_k, with hi_k - lo_k = h and window centres
    spaced by s > h, so at most one section is switching at a time. The bus
    capacitor (C, ESR) integrates the mismatch between delivered section
    current and load current:
        C dVc/dt = i_c,   V_bus = Vc + ESR i_c
    With constant section and load currents every segment is linear, so the
    simulator jumps exactly from one comparator crossing to the next.

Steady state in window j with x = load current not covered by the fully-on
sections (0 < x < Is):
    t_on  = C h_eff / (Is - x),  t_off = C h_eff / x,  h_eff = h - ESR Is
    f     = x (Is - x) / (C h_eff Is),  max at x = Is/2: Is / (4 C h_eff)
    mean V_bus = centre_j + ESR (x - Is/2)
(O'Sullivan and Weinberg, "The sequential switching shunt regulator S3R",
ESA SP-126, 1977, without ESR.)
"""

from dataclasses import dataclass, field
import math

import numpy as np
from scipy.special import lambertw

K_B = 1.380649e-23
Q_E = 1.602176634e-19


# ---------------------------------------------------------------- array model

@dataclass(frozen=True)
class CellSpec:
    """Representative triple-junction GaAs cell at 28 C, AM0, BOL."""
    isc: float = 0.520          # A
    voc: float = 2.700          # V
    n: float = 3.3              # three junctions in series, ~1.1 each
    rs: float = 0.040           # ohm
    t_ref: float = 28.0         # C
    disc_dt: float = 0.32e-3    # A/K
    dvoc_dt: float = -6.0e-3    # V/K
    eol_isc: float = 0.95       # remaining factor after 15 years GEO fluence
    eol_voc: float = 0.94


@dataclass(frozen=True)
class Cell:
    iph: float
    i0: float
    a: float      # n k T / q
    rs: float

    def current(self, v, rs_extra=0.0):
        """Explicit I(V) via Lambert W; rs_extra adds series resistance."""
        rs = self.rs + rs_extra
        v = np.asarray(v, dtype=float)
        s = self.iph + self.i0
        # log of the W argument, to stay finite for large V
        log_arg = math.log(self.i0 * rs / self.a) + (v + s * rs) / self.a
        w = _lambertw_exp(log_arg)
        return s - self.a / rs * w

    def residual(self, v, i, rs_extra=0.0):
        rs = self.rs + rs_extra
        return self.iph - self.i0 * np.expm1((v + i * rs) / self.a) - i


def _lambertw_exp(log_x):
    """W(exp(log_x)), with an asymptotic Newton start where exp overflows."""
    log_x = np.asarray(log_x, dtype=float)
    out = np.empty_like(log_x)
    small = log_x < 500.0
    out[small] = lambertw(np.exp(log_x[small])).real
    big = ~small
    if np.any(big):
        lx = log_x[big]
        w = lx - np.log(lx)
        for _ in range(8):     # solve w + ln w = lx
            w -= (w + np.log(w) - lx) / (1.0 + 1.0 / w)
        out[big] = w
    return out if out.ndim else float(out)


def make_cell(spec: CellSpec, temp_c: float, eol: bool, sun: float = 1.0):
    """Pin Iph and I0 to Isc and Voc at this condition, then scale Iph by sun."""
    isc = (spec.isc + spec.disc_dt * (temp_c - spec.t_ref))
    voc = (spec.voc + spec.dvoc_dt * (temp_c - spec.t_ref))
    if eol:
        isc *= spec.eol_isc
        voc *= spec.eol_voc
    a = spec.n * K_B * (temp_c + 273.15) / Q_E
    # I(0) = Isc and I(Voc) = 0 with Rsh -> infinity:
    #   I0 = Isc / (exp(Voc/a) - exp(Isc Rs/a)),  Iph = I0 (exp(Voc/a) - 1)
    i0 = isc / (math.exp(voc / a) - math.exp(isc * spec.rs / a))
    iph = i0 * math.expm1(voc / a)
    return Cell(iph=iph * sun, i0=i0, a=a, rs=spec.rs)


@dataclass(frozen=True)
class SectionDesign:
    n_series: int = 24
    n_parallel: int = 4
    v_diode: float = 0.8        # blocking diode
    r_harness: float = 0.6      # per string, array to shunt dump, ohm

    def string_current(self, cell: Cell, v_bus):
        v_cell = (np.asarray(v_bus, dtype=float) + self.v_diode) / self.n_series
        return cell.current(v_cell, rs_extra=self.r_harness / self.n_series)

    def section_current(self, cell: Cell, v_bus):
        return self.n_parallel * self.string_current(cell, v_bus)


def max_power_point(cell: Cell):
    v = np.linspace(0.0, 3.2, 32001)
    p = v * cell.current(v)
    k = int(np.argmax(p))
    return v[k], p[k] / v[k]


# ------------------------------------------------------------ regulator model

@dataclass
class S3R:
    i_sec: np.ndarray           # section currents, A (index 0 shunts first)
    c_bus: float                # F
    esr: float                  # ohm
    v0: float                   # centre of window 0, V
    h: float                    # hysteresis width, V
    s: float                    # window spacing, V
    vc: float = 0.0
    shunted: np.ndarray = None
    i_load: float = 0.0
    t: float = 0.0
    log: list = field(default_factory=list)
    trace: list = field(default_factory=list)   # (t, vc, i_cap) per segment

    def __post_init__(self):
        self.i_sec = np.asarray(self.i_sec, dtype=float)
        n = self.i_sec.size
        self.centre = self.v0 + self.s * np.arange(n)
        self.hi = self.centre + self.h / 2
        self.lo = self.centre - self.h / 2
        if self.shunted is None:
            self.shunted = np.zeros(n, dtype=bool)
        if self.esr * self.i_sec.max() >= self.h:
            raise ValueError("ESR step exceeds the hysteresis band: chatter")
        if self.s <= self.h:
            raise ValueError("overlapping windows")

    @property
    def n(self):
        return self.i_sec.size

    def i_cap(self):
        return self.i_sec[~self.shunted].sum() - self.i_load

    def v_bus(self):
        return self.vc + self.esr * self.i_cap()

    def _settle(self):
        """Apply comparator transitions until none is pending."""
        for _ in range(4 * self.n + 4):
            vb = self.v_bus()
            up = np.flatnonzero(~self.shunted & (vb >= self.hi))
            dn = np.flatnonzero(self.shunted & (vb <= self.lo))
            if up.size:
                k = up[0]                 # lowest window crosses first
                self.shunted[k] = True
                self.log.append((self.t, k, +1, vb))
            elif dn.size:
                k = dn[-1]
                self.shunted[k] = False
                self.log.append((self.t, k, -1, vb))
            else:
                return
        raise RuntimeError("comparator cascade did not settle")

    def _next_crossing(self):
        ic = self.i_cap()
        if ic == 0.0:
            return math.inf, None
        slope = ic / self.c_bus
        vb = self.v_bus()
        if slope > 0:
            cand = np.flatnonzero(~self.shunted & (self.hi > vb))
            if not cand.size:
                return math.inf, None
            k = cand[np.argmin(self.hi[cand])]
            return (self.hi[k] - vb) / slope, self.hi[k]
        cand = np.flatnonzero(self.shunted & (self.lo < vb))
        if not cand.size:
            return math.inf, None
        k = cand[np.argmax(self.lo[cand])]
        return (self.lo[k] - vb) / slope, self.lo[k]

    def run(self, t_end, load_steps=()):
        """Advance to t_end. load_steps: sorted (time, current) pairs."""
        steps = [st for st in load_steps if st[0] >= self.t]
        self._settle()
        self.trace.append((self.t, self.vc, self.i_cap()))
        while True:
            dt, target = self._next_crossing()
            t_load = steps[0][0] if steps else math.inf
            t_next = min(self.t + dt, t_load, t_end)
            if t_next == t_end and t_end < min(self.t + dt, t_load):
                self.vc += self.i_cap() / self.c_bus * (t_end - self.t)
                self.t = t_end
                self.trace.append((self.t, self.vc, self.i_cap()))
                return
            if t_load <= self.t + dt:
                self.vc += self.i_cap() / self.c_bus * (t_load - self.t)
                self.t = t_load
                self.i_load = steps.pop(0)[1]
                self.log.append((self.t, -1, 0, self.v_bus()))
            else:
                self.t += dt
                # land exactly on the threshold to avoid round-off misses
                self.vc = target - self.esr * self.i_cap()
            self._settle()
            self.trace.append((self.t, self.vc, self.i_cap()))
            if self.t >= t_end:
                return

    def mean_bus_voltage(self, t1, t2):
        """Exact time average of V_bus over [t1, t2] from the segment trace."""
        tr = self.trace
        total = 0.0
        for (ta, va, ia), (tb, _, _) in zip(tr[:-1], tr[1:]):
            a, b = max(ta, t1), min(tb, t2)
            if b <= a:
                continue
            slope = ia / self.c_bus
            vc_mid = va + slope * ((a + b) / 2 - ta)
            total += (vc_mid + self.esr * ia) * (b - a)
        return total / (t2 - t1)


def steady_state_theory(i_sec, i_load, c, esr, h, s, v0):
    """Active window, deficit, frequency and mean bus voltage (equal sections)."""
    n = len(i_sec)
    i_s = float(i_sec[0])
    full_on = math.floor(i_load / i_s)
    x = i_load - full_on * i_s
    j = n - full_on - 1                  # sections j+1..n-1 fully on
    h_eff = h - esr * i_s
    f = x * (i_s - x) / (c * h_eff * i_s)
    v_mean = v0 + j * s + esr * (x - i_s / 2)
    return dict(j=j, x=x, f=f, duty_on=x / i_s, v_mean=v_mean,
                h_eff=h_eff)


def measure_cycle(sim: S3R, k):
    """Period, on-fraction and mean bus voltage of the last full cycle of k."""
    ev = [(t, a) for t, kk, a, _ in sim.log if kk == k]
    shunt_t = [t for t, a in ev if a == +1]
    rel_t = [t for t, a in ev if a == -1]
    t1, t2 = shunt_t[-2], shunt_t[-1]
    tr = [t for t in rel_t if t1 < t < t2][0]
    period = t2 - t1
    # On-fraction: section k feeds the bus from release to next shunt.
    on_frac = (t2 - tr) / period
    return period, on_frac, sim.mean_bus_voltage(t1, t2)
