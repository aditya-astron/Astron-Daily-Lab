# Astron Built this AP
# Project 046 | Spacecraft power systems | S3R shunt regulator for a 50 V GEO bus
# Author: Aditya (aditya-astron)
# Date: 2026-11-18 IST
"""
Design an S3R for a 50 V regulated GEO bus, then exercise it.

1. String length: the bus must sit on the flat (current-source) part of the
   I-V curve at end of life and hot, so a shorted or released section is a
   clean current step. Rule used: I(V_bus) >= 98 % of Isc at EOL, 70 C.
2. Section count: cover load + battery charge at EOL summer solstice
   (cos 23.44 deg tilt, aphelion 1/1.0167^2) with one section failed.
3. Bus capacitor: worst-case switching frequency is Is/(4 C h_eff) at the
   largest section current (BOL, equinox, hot); size C to keep it <= f_lim.
4. Simulate: frequency sweep over load, the static regulation staircase,
   and a 20 A load step.
"""

import math

import numpy as np

from s3r import (CellSpec, SectionDesign, S3R, make_cell, max_power_point,
                 measure_cycle, steady_state_theory)

V_BUS = 50.0
P_LOAD = 1800.0          # W, platform + payload in sunlight
P_CHARGE = 700.0         # W, battery recharge after the longest eclipse
T_HOT, T_WARM = 70.0, 40.0
SUN_SOLSTICE = math.cos(math.radians(23.44)) / 1.0167 ** 2
SUN_EQUINOX = 1.0
KNEE_RULE = 0.98

H = 0.040                # comparator hysteresis referred to the bus, V
S = 0.050                # window spacing, V
ESR = 0.004              # bus capacitor bank ESR, ohm
F_LIM = 5000.0           # max switching frequency, Hz


def size_string(spec):
    cell = make_cell(spec, T_HOT, eol=True)
    isc = float(cell.current(0.0))
    for ns in range(18, 40):
        sd = SectionDesign(n_series=ns)
        if sd.string_current(cell, V_BUS) >= KNEE_RULE * isc:
            return sd
    raise RuntimeError("no string length meets the knee rule")


def section_current_range(spec, sd):
    cases = {}
    for eol in (False, True):
        for t in (T_WARM, T_HOT):
            for name, sun in (("solstice", SUN_SOLSTICE),
                              ("equinox", SUN_EQUINOX)):
                cell = make_cell(spec, t, eol, sun)
                cases[("EOL" if eol else "BOL", t, name)] = float(
                    sd.section_current(cell, V_BUS))
    return cases


def main():
    spec = CellSpec()
    cell28 = make_cell(spec, 28.0, eol=False)
    vmp, imp = max_power_point(cell28)
    print(f"Cell at 28 C BOL: Vmp = {vmp:.3f} V, Imp = {imp:.4f} A, "
          f"FF = {vmp * imp / (spec.voc * spec.isc):.3f}")

    sd = size_string(spec)
    print(f"\nString: {sd.n_series} cells in series, {sd.n_parallel} strings "
          f"per section (knee rule {KNEE_RULE:.0%} of Isc at EOL {T_HOT:.0f} C)")

    cases = section_current_range(spec, sd)
    for key, i in cases.items():
        print(f"  {key[0]} {key[1]:4.0f} C {key[2]:9s} section current {i:.3f} A")
    i_min, i_max = min(cases.values()), max(cases.values())

    i_req = (P_LOAD + P_CHARGE) / V_BUS
    n_sec = math.ceil(i_req / i_min) + 1
    print(f"\nRequired bus current {i_req:.1f} A, worst section {i_min:.3f} A"
          f" -> {n_sec} sections (N+1)")
    print(f"BOL equinox margin: {n_sec * i_max - i_req:.1f} A shunted at design load"
          f" ({(n_sec * i_max - i_req) / i_max:.1f} sections)")

    # Array I-V slope at the bus: how far from an ideal current source.
    cell_hot_eol = make_cell(spec, T_HOT, eol=True)
    di = (sd.section_current(cell_hot_eol, V_BUS - H / 2)
          - sd.section_current(cell_hot_eol, V_BUS + H / 2))
    print(f"Section current change across the {H * 1e3:.0f} mV band, EOL hot:"
          f" {di * 1e3:.2f} mA ({di / i_min:.1e} of Is)")

    h_eff = H - ESR * i_max
    c_bus = i_max / (4.0 * F_LIM * h_eff)
    c_bus = math.ceil(c_bus / 50e-6) * 50e-6
    f_max = i_max / (4 * c_bus * h_eff)
    print(f"\nBus capacitor: h_eff = {h_eff * 1e3:.2f} mV -> C = {c_bus * 1e3:.2f} mF,"
          f" f_max = {f_max:.0f} Hz at Is = {i_max:.3f} A")
    v0 = V_BUS - S * (n_sec - 1) / 2
    print(f"Regulation band: {v0:.3f} .. {v0 + S * (n_sec - 1):.3f} V,"
          f" sunlight output resistance s/Is = {S / i_max * 1e3:.1f} mohm")

    # Frequency sweep at BOL equinox (largest sections).
    i_sec = np.full(n_sec, i_max)
    print("\nLoad sweep, BOL equinox sections (event-driven simulation vs theory)")
    print("   I_load   win    x/Is   f_sim Hz  f_theory Hz   duty   Vmean sim   Vmean th")
    for frac in (0.08, 0.25, 0.5, 0.75, 0.92):
        i_load = (n_sec // 2 + frac) * i_max
        th = steady_state_theory(i_sec, i_load, c_bus, ESR, H, S, v0)
        sim = S3R(i_sec, c_bus, ESR, v0, H, S, vc=th["v_mean"], i_load=i_load)
        sim.shunted[: th["j"] + 1] = True
        sim.run(20.0 / th["f"])
        period, on_frac, v_mean = measure_cycle(sim, th["j"])
        print(f"  {i_load:7.2f}  {th['j']:4d}  {frac:6.2f}  {1 / period:9.1f}"
              f"  {th['f']:11.1f}  {on_frac:6.3f}  {v_mean:10.5f}  {th['v_mean']:9.5f}")

    # Load step: 20 A on at 2 ms, off at 6 ms.
    i_load0 = 0.35 * n_sec * i_max
    th0 = steady_state_theory(i_sec, i_load0, c_bus, ESR, H, S, v0)
    th1 = steady_state_theory(i_sec, i_load0 + 20.0, c_bus, ESR, H, S, v0)
    sim = S3R(i_sec, c_bus, ESR, v0, H, S, vc=th0["v_mean"], i_load=i_load0)
    sim.shunted[: th0["j"] + 1] = True
    sim.run(10e-3, load_steps=[(2e-3, i_load0 + 20.0), (6e-3, i_load0)])
    vb = np.array([e[3] for e in sim.log])
    tt = np.array([e[0] for e in sim.log])
    on = (tt > 2e-3) & (tt < 6e-3)
    off = tt > 6e-3
    print(f"\n20 A load step from {i_load0:.1f} A: active window {th0['j']} -> {th1['j']}")
    print(f"  min V_bus {vb[on].min():.4f} V (lo of window {th1['j']}:"
          f" {sim.lo[th1['j']]:.4f} V)")
    print(f"  max V_bus after release {vb[off].max():.4f} V (hi of window"
          f" {th0['j']}: {sim.hi[th0['j']]:.4f} V)")
    print(f"  steady deviation {th1['v_mean'] - th0['v_mean']:+.4f} V,"
          f" {sum(1 for e in sim.log if e[1] >= 0)} comparator transitions in 10 ms")


if __name__ == "__main__":
    main()
