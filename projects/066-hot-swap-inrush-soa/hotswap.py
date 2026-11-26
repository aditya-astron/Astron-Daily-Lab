# Astron Built this AP
# Project 066 | Electronics and circuit analysis | Hot-swap inrush control and MOSFET junction temperature
# Author: Aditya (aditya-astron)
# Date: 2026-11-26 IST
"""
dV/dt-controlled hot-swap switch on a 28 V spacecraft bus.

Circuit (high-side N-MOSFET, source follower during the ramp):

    Vb --- D [FET] S --- Vo --+-- C_L --+-- R_L --+-- converter (P, on at V_en)
                G             |         |         |
                +-- Cg -- gnd            gnd       gnd
                +-- Ig (charge pump, clamps at Vb + Vgs_clamp)
                +-- Cgs -- S

Gate and output nodal equations, with the FET as a level-1 square-law device:

    Ig           = Cg dVg + Cgs (dVg - dVo)
    Id(Vgs, Vds) = C_L dVo + Iload(Vo) - Cgs (dVg - dVo)

Once Vgs settles the output slews at S = Ig / Cg and the FET carries
I = C_L S + Iload. For a linear ramp from 0 to Vb (T = Vb / S) the energy
the FET absorbs splits into three closed-form pieces:

    capacitive  E_C = C_L Vb^2 / 2                     (any profile, any S)
    resistive   E_R = Vb^3 / (6 R_L S)
    converter   E_P = (P / S) [Vb ln(Vb / V_en) - (Vb - V_en)]

E_C is exact whenever the output starts at 0 and ends at Vb: the charge
C_L Vb passes through a potential drop that averages Vb / 2. E_R and E_P only
hold for the ideal ramp, so the simulation is the reference for them.

Junction temperature uses a Foster RC network, Zth(t) = sum R_i (1 - e^-t/tau_i).
For power that is piecewise linear in time each Foster state has an exact
one-step update, so Tj needs no extra ODE stiffness.
"""
from dataclasses import dataclass, field, replace

import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class Fet:
    vth: float = 3.0           # V
    k: float = 20.0            # A/V^2, Id = k/2 (Vgs - Vth)^2 in saturation
    cgs: float = 3.0e-9        # F
    # Junction-to-cold-plate Foster network (die, attach, tab, board, plate).
    r_th: tuple = (0.02, 0.15, 0.45, 0.60, 2.80)        # K/W
    tau_th: tuple = (2e-5, 4e-4, 6e-3, 6e-2, 8.0)       # s

    def id(self, vgs, vds):
        vov = vgs - self.vth
        if vov <= 0.0:
            return 0.0
        if vds >= vov:
            return 0.5 * self.k * vov * vov
        return self.k * (vov - 0.5 * vds) * vds

    def rds_on(self, vgs):
        return 1.0 / (self.k * (vgs - self.vth))


@dataclass(frozen=True)
class HotSwap:
    vb: float = 28.0           # bus, V
    c_load: float = 470e-6     # F
    r_load: float = 1.0e3      # bleed / housekeeping, ohm
    p_conv: float = 60.0       # downstream converter input power, W
    v_en: float = 18.0         # converter UVLO turn-on, V
    ig: float = 10e-6          # charge-pump gate current, A
    vgs_clamp: float = 12.0    # gate clamp above the bus, V
    fet: Fet = field(default_factory=Fet)
    pg_gated: bool = False     # converter waits for power-good instead of UVLO

    def cg_for_slew(self, slew):
        """External gate capacitor that sets dVo/dt = slew (V/s)."""
        return self.ig / slew

    def gate_clamp_time(self, cg, ig=None):
        """Time for the gate to reach its clamp (power-good). While the output
        follows, the gate slews at Ig/Cg; while the output is still (below Vth
        at the start, fully on at the end) Cgs loads it too. Vgs sweeps about
        the clamp voltage over those still parts, so
            t = Vb Cg / Ig + Vgs_clamp (Cg + Cgs) / Ig."""
        ig = self.ig if ig is None else ig
        return (self.vb * cg + self.vgs_clamp * (cg + self.fet.cgs)) / ig


# ---------------------------------------------------------------- closed forms

def energy_capacitive(hs):
    return 0.5 * hs.c_load * hs.vb ** 2


def energy_resistive(hs, slew):
    return hs.vb ** 3 / (6.0 * hs.r_load * slew)


def energy_converter(hs, slew):
    if hs.pg_gated:
        return 0.0
    vb, ve = hs.vb, hs.v_en
    return hs.p_conv / slew * (vb * np.log(vb / ve) - (vb - ve))


def energy_turn_on_stall(hs, slew):
    """First-order extra energy in UVLO mode. When the converter starts the FET
    current jumps by P/V_en, so Vgs must rise by
        dVgs = sqrt(2 I+ / k) - sqrt(2 I- / k)
    while the gate keeps slewing at S. The output stalls for ~dVgs / S with
    (Vb - V_en) across the FET and the converter current through it."""
    if hs.pg_gated:
        return 0.0
    k = hs.fet.k
    i_minus = hs.c_load * slew + hs.v_en / hs.r_load
    i_plus = i_minus + hs.p_conv / hs.v_en
    dvgs = np.sqrt(2 * i_plus / k) - np.sqrt(2 * i_minus / k)
    return (hs.vb - hs.v_en) * hs.p_conv / hs.v_en * dvgs / slew


def energy_closed_form(hs, slew):
    return energy_capacitive(hs) + energy_resistive(hs, slew) + energy_converter(hs, slew)


def peak_current_closed_form(hs, slew):
    """Largest FET current on the ideal ramp: C_L S + Vo/R_L + P/Vo, worst at
    the converter turn-on point or at the top of the ramp."""
    base = hs.c_load * slew
    if hs.pg_gated:
        # ramp current, or the full load stepping on after power-good
        return max(base, hs.p_conv / hs.vb) + hs.vb / hs.r_load
    return base + max(hs.v_en / hs.r_load + hs.p_conv / hs.v_en,
                      hs.vb / hs.r_load + hs.p_conv / hs.vb)


# ------------------------------------------------------------------ simulation

@dataclass
class Trace:
    t: np.ndarray
    vg: np.ndarray
    vo: np.ndarray
    id: np.ndarray
    p: np.ndarray            # FET dissipation, W
    t_enable: float          # converter turn-on time
    t_done: float            # gate reached its clamp

    @property
    def energy(self):
        return np.trapezoid(self.p, self.t)


def simulate(hs, slew, n_out=20001, settle=5e-3):
    """Integrate the ramp. Phase 1: converter off. Phase 2: converter on
    (latched; its UVLO hysteresis is wider than any dip seen here) until the
    gate clamps. Phase 3: gate clamped, output settles under full load."""
    fet = hs.fet
    cg = hs.cg_for_slew(slew)
    cgs = fet.cgs
    vclamp = hs.vb + hs.vgs_clamp
    m = np.array([[cg + cgs, -cgs], [-cgs, hs.c_load + cgs]])
    m_inv = np.linalg.inv(m)

    def iload(vo, conv_on):
        i = vo / hs.r_load
        if conv_on:
            i += hs.p_conv / max(vo, 1e-3)
        return i

    def rhs(_t, y, conv_on, clamped):
        vg, vo = y
        idr = fet.id(vg - vo, hs.vb - vo)
        il = iload(vo, conv_on)
        if clamped:
            return [0.0, (idr - il) / (hs.c_load + cgs)]
        return m_inv @ np.array([hs.ig, idr - il])

    def ev_enable(_t, y, *_):
        return (y[0] - (vclamp - 1e-6)) if hs.pg_gated else (y[1] - hs.v_en)
    ev_enable.terminal = True
    ev_enable.direction = 1

    def ev_clamp(_t, y, *_):
        return y[0] - vclamp
    ev_clamp.terminal = True
    ev_clamp.direction = 1

    t_max = 50.0 * hs.vb / slew + 1.0
    opts = dict(method="Radau", rtol=1e-9, atol=1e-12, dense_output=True)
    segs = []

    s1 = solve_ivp(rhs, (0.0, t_max), [0.0, 0.0], args=(False, False),
                   events=ev_enable, **opts)
    t_en = s1.t[-1]
    segs.append((s1, False, False))
    if hs.pg_gated:
        t_done = t_en
        y_cl = [vclamp, s1.y[1, -1]]
    else:
        s2 = solve_ivp(rhs, (t_en, t_max), s1.y[:, -1], args=(True, False),
                       events=ev_clamp, **opts)
        segs.append((s2, True, False))
        t_done = s2.t[-1]
        y_cl = [vclamp, s2.y[1, -1]]
    s3 = solve_ivp(rhs, (t_done, t_done + settle), y_cl, args=(True, True), **opts)
    segs.append((s3, True, True))

    # Uniform sampling inside each segment plus both edges, so the FET power
    # step at converter turn-on is kept and integrated exactly.
    total = segs[-1][0].t[-1]
    ts, vgs, vos = [], [], []
    for sol, _, _ in segs:
        a, b = sol.t[0], sol.t[-1]
        n = max(int(n_out * (b - a) / total), 200)
        tt = np.linspace(a, b, n)
        yy = sol.sol(tt)
        ts.append(tt)
        vgs.append(yy[0])
        vos.append(yy[1])
    t = np.concatenate(ts)
    vg = np.concatenate(vgs)
    vo = np.concatenate(vos)
    idr = np.array([fet.id(g - o, hs.vb - o) for g, o in zip(vg, vo)])
    p = idr * (hs.vb - vo)
    return Trace(t, vg, vo, idr, p, t_en, t_done)


# ------------------------------------------------------------- thermal network

def zth(fet, t):
    t = np.asarray(t, dtype=float)
    r = np.asarray(fet.r_th)[:, None]
    tau = np.asarray(fet.tau_th)[:, None]
    return np.sum(r * (1.0 - np.exp(-np.atleast_1d(t)[None, :] / tau)), axis=0)


def junction_rise(fet, t, p, t_cool=0.2, n_cool=2000):
    """Foster-network temperature rise for power p(t) linear between samples,
    followed by t_cool seconds with p = 0. Exact per step:
        theta' = (R P(t) - theta) / tau,  P = P0 + m s on [0, h]
        theta(h) = theta0 a + R (P1 - P0 a) - R (P1 - P0) (1 - a) tau / h,
        a = e^-h/tau, written with expm1 so a near-zero step stays exact.
    Returns (time, rise)."""
    if t_cool > 0:
        tc = t[-1] + np.linspace(0.0, t_cool, n_cool + 1)[1:]
        # one sample just after the end so power drops to zero in ~0 time
        t = np.concatenate([t, [t[-1] + 1e-9], tc])
        p = np.concatenate([p, [0.0], np.zeros(n_cool)])
    r = np.asarray(fet.r_th)
    tau = np.asarray(fet.tau_th)
    theta = np.zeros_like(r)
    rise = np.zeros_like(t)
    for k in range(1, len(t)):
        h = t[k] - t[k - 1]
        if h <= 0:
            rise[k] = rise[k - 1]
            continue
        p0, p1 = p[k - 1], p[k]
        x = h / tau
        a = np.exp(-x)
        g = -np.expm1(-x) / x          # (1 - a) tau / h without cancellation
        theta = theta * a + r * (p1 - p0 * a) - r * (p1 - p0) * g
        rise[k] = theta.sum()
    return t, rise


def peak_tj(hs, slew, t_amb=60.0, **kw):
    tr = simulate(hs, slew, **kw)
    _, rise = junction_rise(hs.fet, tr.t, tr.p)
    return t_amb + rise.max(), tr


def with_(hs, **kw):
    return replace(hs, **kw)
