# Astron Built this AP
# Project 050 | Rocket propulsion | Regenerative cooling of a LOX/methane thrust chamber
# Author: Aditya (aditya-astron)
# Date: 2026-11-20 IST
"""
Design case: 7 MPa LOX/methane chamber, 60 mm throat, all fuel through a
milled CuCrZr jacket in counterflow from the nozzle exit.

1. Baseline jacket: wall temperature, heat flux and coolant state along the
   chamber, total heat load and pressure drop.
2. 2D unit-cell conduction check of the 1D wall/fin model at the throat.
3. Trade over channel count and tapered channel height against a 900 K
   hot-wall limit and a jacket pressure-drop budget of 25 % of chamber
   pressure.
"""
from __future__ import annotations

import math

import numpy as np

from regen import (ChamberContour, CombustionGas, CoolPropFluid, Jacket, RegenSolver,
                   bartz_gas_side, solve_station)
from unitcell import fin_model_flat, solve_cell

# Combustion products, LOX/CH4 at O/F ~ 3.4 and 7 MPa (frozen chamber values of
# the order given by an equilibrium code; treated as inputs here).
GAS = CombustionGas(T0=3550.0, gamma=1.15, molar_mass=21.0, pc=7.0e6)
OF = 3.4
CONTOUR = ChamberContour(Rt=0.030, contraction=4.0, lstar=1.0, eps_exit=8.0)
T_IN, P_IN = 120.0, 12.0e6         # coolant at the jacket inlet manifold
TWG_LIMIT = 900.0                  # CuCrZr hot-wall design limit, K
DP_LIMIT = 0.25 * GAS.pc


def taper(h_throat: float, h_chamber: float):
    """Channel height linear in radius from the throat to the chamber radius,
    held at h_chamber beyond it (barrel and the wide end of the nozzle)."""
    Rt, Rc = CONTOUR.Rt, CONTOUR.Rc
    return lambda r: h_throat + (h_chamber - h_throat) * min(max((r - Rt) / (Rc - Rt), 0.0), 1.0)


def baseline_jacket(n=120, h_throat=3.0e-3, h_chamber=1.5e-3):
    return Jacket(n_channels=n, height=taper(h_throat, h_chamber), land=0.8e-3,
                  t_wall=0.7e-3, k_wall=320.0)


def engine_flows():
    At = math.pi * CONTOUR.Rt ** 2
    mdot = GAS.pc * At / GAS.cstar
    return mdot, mdot / (1.0 + OF)


def main():
    fluid = CoolPropFluid("Methane")
    gas_side = bartz_gas_side(GAS, CONTOUR)
    mdot, mdot_f = engine_flows()

    print("Engine")
    print(f"  c* {GAS.cstar:7.1f} m/s   cp {GAS.cp:6.1f} J/kgK   Pr {GAS.prandtl:.3f}"
          f"   mu0 {GAS.mu0 * 1e5:.3f}e-5 Pa s")
    print(f"  mdot {mdot:.3f} kg/s, fuel {mdot_f:.3f} kg/s, Rt {CONTOUR.Rt * 1e3:.0f} mm,"
          f" Rc {CONTOUR.Rc * 1e3:.0f} mm, Re {CONTOUR.Re * 1e3:.1f} mm")
    print(f"  cylinder {CONTOUR.Lcyl * 1e3:.1f} mm, injector at x = {CONTOUR.x_inj * 1e3:.1f} mm,"
          f" exit at x = {CONTOUR.x_exit * 1e3:.1f} mm")

    jk = baseline_jacket()
    sol = RegenSolver(CONTOUR, jk, fluid, mdot_f, gas_side).run(T_IN, P_IN, n=400)
    i_max = int(np.argmax(sol.Twg))
    i_q = int(np.argmax(sol.q))
    print(f"\nBaseline jacket: N = {jk.n_channels}, h = 3.0 mm at throat -> 1.5 mm at chamber,"
          f" land {jk.land * 1e3:.1f} mm, ligament {jk.t_wall * 1e3:.1f} mm")
    print("     x[mm]  r[mm]     M   hg[kW/m2K]  q[MW/m2]  Twg[K]  Twc[K]  Tb[K]  p[MPa]")
    for xm in (CONTOUR.x_exit, 0.10, 0.05, 0.0, -0.02, -0.06, CONTOUR.x_inj):
        j = int(np.argmin(abs(sol.x - xm)))
        print(f"  {sol.x[j] * 1e3:8.1f} {sol.r[j] * 1e3:6.1f} {sol.M[j]:5.2f} {sol.hg[j] / 1e3:10.2f}"
              f" {sol.q[j] / 1e6:9.2f} {sol.Twg[j]:7.1f} {sol.Twc[j]:7.1f} {sol.Tb[j]:6.1f}"
              f" {sol.p[j] / 1e6:7.3f}")
    print(f"  peak q {sol.q[i_q] / 1e6:.2f} MW/m2 at x = {sol.x[i_q] * 1e3:.1f} mm")
    print(f"  peak Twg {sol.Twg[i_max]:.1f} K at x = {sol.x[i_max] * 1e3:.1f} mm")
    print(f"  heat load {sol.Q_total / 1e6:.3f} MW, coolant {sol.Tb[0]:.0f} -> {sol.Tb[-1]:.1f} K,"
          f" dp {sol.dp / 1e6:.3f} MPa")

    # --- 2D check at the throat -------------------------------------------------
    j = int(np.argmin(abs(sol.x)))
    x0, r0 = sol.x[j], sol.r[j]
    st = fluid.state(sol.h[j], sol.p[j])
    res = solve_station(jk, r0, lambda T: sol.hg[j], gas_side(x0, 1000.0)[1], st,
                        mdot_f / jk.n_channels, lambda T: fluid.viscosity(T, sol.p[j]))
    w, b0, h0, _, _ = jk.channel(r0)
    cell = (w, b0, jk.t_wall, h0, jk.k_wall, sol.hg[j], res.Taw, res.hc, st.T)
    Q1, T1 = fin_model_flat(*cell)
    print("\n2D unit cell at the throat vs. the 1D fin model on the same flat cell")
    print(f"  w {w * 1e3:.3f} mm, land {b0 * 1e3:.3f} mm, h {h0 * 1e3:.2f} mm,"
          f" hg {sol.hg[j] / 1e3:.2f} kW/m2K, hc {res.hc / 1e3:.2f} kW/m2K, Tb {st.T:.1f} K")
    print(f"  1D: Twg {T1:.1f} K (annular march gave {res.Twg:.1f} K)")
    print("    n   Q2D/Q1D   Twg max 2D   Twg min 2D")
    for n in (10, 20, 40, 80):
        Q2, Ts, _, _ = solve_cell(*cell, n=n)
        print(f"  {n:3d}   {Q2 / Q1:7.4f}   {Ts.max():10.1f}   {Ts.min():10.1f}")

    # --- trade ------------------------------------------------------------------
    print(f"\nTrade (limits: Twg <= {TWG_LIMIT:.0f} K, dp <= {DP_LIMIT / 1e6:.1f} MPa)")
    print("    N  h_t/h_c[mm]  Twg max[K] at x[mm]  Twg throat[K]  dp[MPa]  Tout[K]  ok")
    best = None
    for n in (96, 108, 120, 132):
        for ht, hc in ((2.5e-3, 1.5e-3), (3.0e-3, 1.5e-3), (2.5e-3, 1.8e-3), (3.0e-3, 1.8e-3)):
            tag = f"{ht * 1e3:.1f}/{hc * 1e3:.1f}"
            try:
                s = RegenSolver(CONTOUR, baseline_jacket(n, ht, hc), fluid, mdot_f,
                                gas_side).run(T_IN, P_IN, n=200)
            except ValueError as e:
                print(f"  {n:3d}  {tag:>11}   infeasible ({e})")
                continue
            i = int(np.argmax(s.Twg))
            ok = s.Twg[i] <= TWG_LIMIT and s.dp <= DP_LIMIT
            print(f"  {n:3d}  {tag:>11}  {s.Twg[i]:10.1f} {s.x[i] * 1e3:7.1f}"
                  f"  {s.Twg[np.argmin(abs(s.x))]:13.1f}  {s.dp / 1e6:7.3f}  {s.Tb[-1]:7.1f}"
                  f"  {'yes' if ok else 'no'}")
            if ok and (best is None or s.dp < best[2]):
                best = (n, tag, s.dp, s.Twg[i])
    if best:
        print(f"  lowest-dp feasible: N = {best[0]}, h = {best[1]} mm,"
              f" dp {best[2] / 1e6:.3f} MPa, Twg max {best[3]:.1f} K")
    else:
        print("  no feasible design in the grid")


if __name__ == "__main__":
    main()
