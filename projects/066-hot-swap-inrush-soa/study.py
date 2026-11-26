# Astron Built this AP
# Project 066 | Electronics and circuit analysis | Hot-swap inrush control and MOSFET junction temperature
# Author: Aditya (aditya-astron)
# Date: 2026-11-26 IST
"""
Size the gate capacitor of a dV/dt-controlled hot-swap switch that feeds a
120 W payload converter (470 uF input bank) from a 28 V bus, with the FET
on a 60 C cold plate and junction derated to 110 C. Compare letting the
converter start at its 18 V UVLO against holding it off until power-good.
"""
import numpy as np
from scipy.optimize import brentq, minimize_scalar

import hotswap as h

T_AMB = 60.0
TJ_MAX = 110.0
P_PAYLOAD = 120.0
T_BUDGET = 0.100           # gate must clamp (power-good) within 100 ms
IG_TOL = 0.20              # charge-pump current tolerance
CG_TOL = 0.10              # capacitor tolerance
E6 = [4.7e-9, 6.8e-9, 10e-9, 15e-9, 22e-9, 33e-9, 47e-9]


def hand_check(hs):
    s = 1000.0
    print("Hand check, S = 1 V/ms, capacitive load only")
    print(f"  Cg = Ig/S = {hs.cg_for_slew(s) * 1e9:.1f} nF, ramp Vb/S = {hs.vb / s * 1e3:.1f} ms")
    print(f"  I = C S = {hs.c_load * s:.3f} A, peak FET power C S Vb = {hs.c_load * s * hs.vb:.2f} W")
    print(f"  E = C Vb^2/2 = {h.energy_capacitive(hs):.5f} J")
    x = h.with_(hs, p_conv=0.0, r_load=1e12)
    tr = h.simulate(x, s)
    print(f"  simulated  {tr.energy:.5f} J, final Vo {tr.vo[-1]:.6f} V\n")


def energy_table(hs):
    print("FET energy: simulation vs ramp closed forms (J)")
    print("  case                       S V/s    sim       E_C+E_R+E_P  +stall    sim/cf")
    cases = [
        ("C only", dict(p_conv=0.0, r_load=1e12)),
        ("C + 50 ohm", dict(p_conv=0.0, r_load=50.0)),
        ("C + 120 W, power-good", dict(pg_gated=True)),
        ("C + 120 W, UVLO 18 V", dict()),
    ]
    for name, kw in cases:
        x = h.with_(hs, **kw)
        for s in (300.0, 3000.0):
            tr = h.simulate(x, s)
            cf = h.energy_closed_form(x, s)
            st = h.energy_turn_on_stall(x, s)
            print(f"  {name:24s} {s:7.0f}  {tr.energy:8.5f}  {cf:10.5f}  {cf + st:8.5f}  {tr.energy / cf:7.4f}")
    print()


def thermal_check(fet):
    print("Foster network check")
    for tp in (1e-4, 1e-3, 1e-2, 1e-1):
        t = np.array([0.0, 1e-12, tp])
        p = np.array([0.0, 10.0, 10.0])
        _, rise = h.junction_rise(fet, t, p, t_cool=0)
        print(f"  10 W for {tp * 1e3:6.1f} ms: stepped {rise[-1]:.6f} K, P Zth(t) {10 * h.zth(fet, tp)[0]:.6f} K")
    print(f"  Zth(inf) = {sum(fet.r_th):.2f} K/W junction to plate\n")


def sweep(hs):
    print(f"Slew sweep, {P_PAYLOAD:.0f} W converter, plate {T_AMB:.0f} C")
    print("    S V/s   Cg nF   t_pg ms |  UVLO: E J   Ipk A   Tj C   dip V |  PG: E J   Ipk A   Tj C")
    for s in (100, 200, 500, 1000, 2000, 5000, 10000, 20000):
        out = []
        for pg in (False, True):
            x = h.with_(hs, pg_gated=pg)
            tj, tr = h.peak_tj(x, s, T_AMB, n_out=4001)
            dip = max(hs.v_en - tr.vo[tr.t >= tr.t_enable].min(), 0.0) if not pg else 0.0
            out.append((tr.energy, tr.id.max(), tj, tr.t_done, dip))
        print(f"  {s:7.0f}  {hs.cg_for_slew(s) * 1e9:6.2f}  {out[1][3] * 1e3:8.1f} |"
              f"  {out[0][0]:7.3f}  {out[0][1]:6.2f}  {out[0][2]:6.1f}  {out[0][4]:6.2f} |"
              f"  {out[1][0]:7.3f}  {out[1][1]:6.2f}  {out[1][2]:6.1f}")
    print()


def best_slew_uvlo(hs):
    f = lambda ls: h.peak_tj(hs, 10 ** ls, T_AMB, n_out=3001)[0]
    r = minimize_scalar(f, bounds=(2.0, 4.7), method="bounded", options={"xatol": 0.01})
    return 10 ** r.x, r.fun


def power_limits(hs):
    print("Largest converter that starts within Tj 110 C")
    for p in (60.0, 120.0, 200.0):
        s, tj = best_slew_uvlo(h.with_(hs, p_conv=p))
        print(f"  UVLO start, P = {p:5.0f} W: best S = {s:6.0f} V/s, Tj = {tj:6.1f} C")

    def margin(p):
        return best_slew_uvlo(h.with_(hs, p_conv=p))[1] - TJ_MAX
    p_max = brentq(margin, 120.0, 400.0, xtol=0.5)
    x = h.with_(hs, p_conv=p_max)
    s_opt, _ = best_slew_uvlo(x)
    ipk = h.simulate(x, s_opt).id.max()
    print(f"  UVLO start limit: {p_max:.0f} W (at S = {s_opt:.0f} V/s, inrush peak {ipk:.1f} A)")
    x = h.with_(hs, pg_gated=True, p_conv=1000.0)
    tj, _ = h.peak_tj(x, 1000.0, T_AMB, n_out=4001)
    print(f"  Power-good start: ramp heating does not depend on P; even 1 kW at 1 V/ms gives Tj = {tj:.1f} C\n")


def pick_cg(hs):
    print(f"Gate capacitor pick, power-good gated, gate clamp within {T_BUDGET * 1e3:.0f} ms at the slow corner")
    chosen = None
    for cg in E6:
        s_slow = hs.ig * (1 - IG_TOL) / (cg * (1 + CG_TOL))
        s_fast = hs.ig * (1 + IG_TOL) / (cg * (1 - CG_TOL))
        t_slow = hs.gate_clamp_time(cg * (1 + CG_TOL), hs.ig * (1 - IG_TOL))
        ok = t_slow <= T_BUDGET
        print(f"  Cg {cg * 1e9:5.1f} nF: S {s_slow:6.0f}..{s_fast:6.0f} V/s, gate clamp at slow corner {t_slow * 1e3:6.1f} ms {'ok' if ok else 'late'}")
        if ok:
            chosen = (cg, s_slow, s_fast)
    cg, s_slow, s_fast = chosen
    pg = h.with_(hs, pg_gated=True)
    tj_fast, trf = h.peak_tj(pg, s_fast, T_AMB)
    tj_uv, tru = h.peak_tj(hs, s_fast, T_AMB)
    print(f"  Chosen Cg = {cg * 1e9:.0f} nF")
    print(f"  fast corner S = {s_fast:.0f} V/s: power-good Tj {tj_fast:.1f} C, Ipk {trf.id.max():.2f} A, E {trf.energy:.3f} J")
    print(f"                       UVLO start Tj {tj_uv:.1f} C, Ipk {tru.id.max():.2f} A, E {tru.energy:.3f} J")
    tj_slow_uv, _ = h.peak_tj(hs, s_slow, T_AMB)
    print(f"  slow corner S = {s_slow:.0f} V/s: UVLO start Tj {tj_slow_uv:.1f} C")


def main():
    hs = h.HotSwap(p_conv=P_PAYLOAD)
    hand_check(hs)
    energy_table(hs)
    thermal_check(hs.fet)
    sweep(hs)
    power_limits(hs)
    pick_cg(hs)


if __name__ == "__main__":
    main()
