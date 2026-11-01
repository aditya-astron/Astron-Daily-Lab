# Astron Built this AP
# Project 004 | Spacecraft power systems | Solar array MPPT through a LEO eclipse cycle
# Author: Aditya (aditya-astron)
# Date: 2026-11-01 IST
"""Simulate a sun-pointed solar array through repeated 500 km LEO orbits and
compare three ways of loading it:

  DET     direct energy transfer, array clamped to a fixed bus voltage
  P&O     perturb and observe MPPT
  IncCond incremental conductance MPPT

The interesting part is eclipse exit. The panel comes out of shadow near
-80 C, where Vmp is ~30 % above its hot value. A fixed-voltage array wastes
that margin; a tracker has to walk to it and then follow it back down as the
panel warms over the next ~20 minutes.
"""

import math
from dataclasses import dataclass

import numpy as np

from pv_model import (G_AM0, Array, CellDatasheet, extract_params, params_at,
                      true_mpp)

MU_EARTH = 398600.4418e9   # m^3/s^2
R_EARTH = 6378.137e3       # m
SIGMA = 5.670374419e-8     # W/m^2/K^4


@dataclass(frozen=True)
class Scenario:
    alt_m: float = 500e3
    beta_deg: float = 0.0          # sun angle to orbit plane; 0 = longest eclipse
    penumbra_s: float = 9.0        # linear irradiance ramp at shadow edges
    g_sun: float = G_AM0           # W/m^2
    dt: float = 0.1                # s, also the MPPT update period
    n_orbits: int = 3
    # Panel thermal model (per m^2 of panel): single node, radiating from
    # both faces, Earth IR absorbed on the back. Albedo is ignored.
    absorptance: float = 0.91
    emit_front: float = 0.85
    emit_back: float = 0.80
    q_earth_ir: float = 70.0       # W/m^2 absorbed
    heat_cap: float = 2500.0       # J/m^2/K, ~2.7 kg/m^2 honeycomb panel
    v_det: float = 35.0            # V, fixed array voltage for DET
    v_park: float = 40.0           # V, where trackers sit during eclipse
    mppt_step: float = 0.2         # V per update for both trackers


def orbit_period(alt_m):
    a = R_EARTH + alt_m
    return 2.0 * math.pi * math.sqrt(a ** 3 / MU_EARTH)


def eclipse_half_angle(alt_m, beta_deg):
    """Cylindrical shadow. The satellite is in shadow when the in-plane
    angle from the anti-sun direction is below this value."""
    r = R_EARTH + alt_m
    beta = math.radians(beta_deg)
    c = math.sqrt(max(0.0, 1.0 - (R_EARTH / r) ** 2)) / math.cos(beta)
    if c >= 1.0:
        return 0.0                 # no eclipse at this beta
    return math.acos(c)


def irradiance_profile(sc: Scenario):
    """Irradiance on the panel vs time. t = 0 is eclipse exit, so every orbit
    starts with the cold-panel transient."""
    period = orbit_period(sc.alt_m)
    half = eclipse_half_angle(sc.alt_m, sc.beta_deg)
    t_ecl = period * (2.0 * half) / (2.0 * math.pi)
    t = np.arange(0.0, sc.n_orbits * period, sc.dt)
    phase = np.mod(t, period)
    t_sun = period - t_ecl
    # Signed time from the nearer shadow edge (positive in sunlight), mapped to
    # a ramp centred on the umbra boundary of the cylindrical model.
    in_sun = phase <= t_sun
    edge = np.where(in_sun, np.minimum(phase, t_sun - phase),
                    -np.minimum(phase - t_sun, period - phase))
    g = sc.g_sun * np.clip(edge / sc.penumbra_s + 0.5, 0.0, 1.0)
    return t, g, period, t_ecl


def panel_temperature(g, sc: Scenario, eff_elec=0.29):
    """Single-node energy balance per m^2:
        C dT/dt = alpha G + q_IR - (eps_f + eps_b) sigma T^4 - eta G
    Electrical power leaves the panel, so it is subtracted using the cell
    efficiency at MPP. Integrated with forward Euler (tau ~ 10 min >> dt).
    The start temperature is iterated until the orbit is periodic."""
    eps = sc.emit_front + sc.emit_back
    q_in = (sc.absorptance - eff_elec) * g + sc.q_earth_ir
    t_k = 250.0
    for _ in range(6):
        temps = np.empty_like(g)
        tk = t_k
        for i, q in enumerate(q_in):
            temps[i] = tk
            tk += sc.dt * (q - eps * SIGMA * tk ** 4) / sc.heat_cap
        if abs(tk - t_k) < 1e-3:
            break
        t_k = tk
    return temps - 273.15


class PerturbObserve:
    """Fixed-step hill climb on measured power."""

    def __init__(self, v0, step=0.2):
        self.v = v0
        self.step = step
        self.p_prev = None
        self.v_prev = v0
        self.direction = 1.0

    def update(self, v, i):
        p = v * i
        if self.p_prev is None:
            self.p_prev, self.v_prev = p, v
            self.v = v + self.direction * self.step
            return self.v
        dp, dv = p - self.p_prev, v - self.v_prev
        if dp != 0.0:
            # Keep going if the last move raised power, otherwise reverse.
            self.direction = math.copysign(1.0, dp) * math.copysign(1.0, dv)
        self.p_prev, self.v_prev = p, v
        self.v = v + self.direction * self.step
        return self.v


class IncrementalConductance:
    """At the MPP dP/dV = I + V dI/dV = 0, i.e. dI/dV = -I/V. The sign of
    dI/dV + I/V says which side of the peak we are on, without needing the
    previous move to have been in a known direction."""

    def __init__(self, v0, step=0.2, tol=0.005):
        self.v = v0
        self.step = step
        self.tol = tol
        self.v_prev = v0
        self.i_prev = None

    def update(self, v, i):
        if self.i_prev is None:
            # No previous sample yet (start, or just out of eclipse): a
            # difference against the dark reading would be all irradiance
            # change and no V change. Take one plain step to get a baseline.
            self.v_prev, self.i_prev = v, i
            self.v = v + self.step
            return self.v
        dv, di = v - self.v_prev, i - self.i_prev
        self.v_prev, self.i_prev = v, i
        if v <= 0.0 or i <= 0.0:
            move = -1.0                  # past Voc (or no light): come down
        elif dv == 0.0:
            # Holding still: a current change means irradiance moved the MPP.
            move = 0.0 if di == 0.0 else math.copysign(1.0, di)
        else:
            s = di / dv + i / v          # proportional to dP/dV
            move = 0.0 if abs(s) < self.tol * i / v else math.copysign(1.0, s)
        self.v = v + move * self.step
        return self.v


def run(sc: Scenario = Scenario(), noise=True, seed=4):
    ds = CellDatasheet()
    ref = extract_params(ds)
    array = Array()

    t, g, period, t_ecl = irradiance_profile(sc)
    temp_c = panel_temperature(g, sc)
    cond = [params_at(ref, ds, gi, ti) for gi, ti in zip(g, temp_c)]
    v_mpp, p_mpp = true_mpp(array, cond)

    rng = np.random.default_rng(seed)
    # Sensor noise, 1 sigma: 20 mV and 10 mA (roughly 12-bit ADCs on a
    # 60 V / 25 A range after filtering).
    sv, si = (0.02, 0.01) if noise else (0.0, 0.0)

    def det_power(c):
        return float(array.power(sc.v_det, c)) if c.iph > 0 else 0.0

    trackers = {
        "P&O": PerturbObserve(sc.v_park, sc.mppt_step),
        "IncCond": IncrementalConductance(sc.v_park, sc.mppt_step),
    }
    p_out = {"DET": np.array([max(det_power(c), 0.0) for c in cond])}
    v_log = {}
    for name, trk in trackers.items():
        v_op = sc.v_park
        pw = np.empty_like(g)
        vl = np.empty_like(g)
        for k, c in enumerate(cond):
            i_true = float(array.current(v_op, c)) if c.iph > 0 else 0.0
            i_true = max(i_true, 0.0)     # blocking diode stops reverse flow
            pw[k] = v_op * i_true
            vl[k] = v_op
            if g[k] <= 0.0:
                # Eclipse: no information in the measurements, so park.
                v_op = sc.v_park
                trk.__init__(sc.v_park, sc.mppt_step)
                continue
            v_meas = v_op + rng.normal(0.0, sv)
            i_meas = i_true + rng.normal(0.0, si)
            v_cmd = trk.update(v_meas, i_meas)
            v_op = min(max(v_cmd, 0.0), 60.0)
        p_out[name] = pw
        v_log[name] = vl

    return dict(t=t, g=g, temp_c=temp_c, v_mpp=v_mpp, p_mpp=p_mpp,
                p_out=p_out, v_log=v_log, period=period, t_ecl=t_ecl, sc=sc)


def settle_time(t, p, p_ref, start, frac=0.99, hold=5.0):
    """Seconds after `start` until p stays above frac * p_ref for `hold` s.
    Only lit samples count; in eclipse 0 >= 0.99 * 0 would pass trivially."""
    mask = (t >= start) & (p_ref > 0)
    tt, ok = t[mask], (p[mask] >= frac * p_ref[mask])
    n_hold = int(round(hold / (tt[1] - tt[0])))
    run_len = 0
    for k, good in enumerate(ok):
        run_len = run_len + 1 if good else 0
        if run_len >= n_hold:
            return tt[k - n_hold + 1] - start
    return math.inf


def summarize(res):
    t, period, dt = res["t"], res["period"], res["sc"].dt
    last = t >= (res["sc"].n_orbits - 1) * period   # periodic steady state
    e_avail = res["p_mpp"][last].sum() * dt / 3600.0
    sunlit = last & (res["g"] >= res["sc"].g_sun)
    print(f"Orbit: {res['sc'].alt_m/1e3:.0f} km, beta {res['sc'].beta_deg:.0f} deg, "
          f"period {period/60:.2f} min, eclipse {res['t_ecl']/60:.2f} min")
    print(f"Panel temperature over orbit: {res['temp_c'][last].min():.1f} C "
          f"to {res['temp_c'][last].max():.1f} C")
    print(f"Vmp in full sun: {res['v_mpp'][sunlit].min():.2f} V "
          f"to {res['v_mpp'][sunlit].max():.2f} V")
    print(f"Pmp in full sun: {res['p_mpp'][sunlit].min():.1f} W "
          f"to {res['p_mpp'][sunlit].max():.1f} W")
    print(f"Energy available at MPP: {e_avail:.2f} Wh per orbit\n")
    print(f"{'method':<9}{'Wh/orbit':>10}{'efficiency':>12}{'99% settle':>12}")
    start = (res["sc"].n_orbits - 1) * period
    out = {}
    for name, p in res["p_out"].items():
        e = p[last].sum() * dt / 3600.0
        ts = settle_time(t, p, res["p_mpp"], start)
        ts_txt = "never" if math.isinf(ts) else f"{ts:.1f} s"
        print(f"{name:<9}{e:>10.2f}{100*e/e_avail:>11.2f}%{ts_txt:>12}")
        out[name] = e / e_avail
    return out


def step_sweep(steps=(0.02, 0.05, 0.2, 0.5, 1.0, 2.0)):
    """P&O step size trade: small steps are slow to reach the cold Vmp at
    eclipse exit, large steps waste power oscillating around the peak."""
    print(f"\nP&O step size sweep (10 Hz updates)")
    print(f"{'step V':>7}{'efficiency':>12}{'99% settle':>12}")
    rows = []
    for s in steps:
        res = run(Scenario(mppt_step=s))
        period, dt, t = res["period"], res["sc"].dt, res["t"]
        last = t >= (res["sc"].n_orbits - 1) * period
        eff = res["p_out"]["P&O"][last].sum() / res["p_mpp"][last].sum()
        ts = settle_time(t, res["p_out"]["P&O"], res["p_mpp"],
                         (res["sc"].n_orbits - 1) * period)
        ts_txt = "never" if math.isinf(ts) else f"{ts:.1f} s"
        print(f"{s:>7.2f}{100*eff:>11.3f}%{ts_txt:>12}")
        rows.append((s, eff, ts))
    return rows


if __name__ == "__main__":
    summarize(run())
    step_sweep()
