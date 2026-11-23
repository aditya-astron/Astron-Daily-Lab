# Astron Built this AP
# Project 060 | Navigation / Guidance | Lunar terminal descent with energy-optimal ZEM/ZEV guidance
# Author: Aditya (aditya-astron)
# Date: 2026-11-23 IST
"""Run the descent study: time-to-go trade, a nominal landing, and limits."""
import numpy as np
from landing import (G0, G_MOON, Lander, energy, feasible_tgo_window, ideal_dv,
                     ideal_peak_accel, optimal_tgo, simulate, profile_coeffs)

R0 = np.array([-2000.0, 400.0, 1500.0])   # m, target-centred (x downrange, z up)
V0 = np.array([45.0, -5.0, -25.0])        # m/s


def summary(label, res, lander):
    m = res.x[:, 6]
    dv = res.x[-1, 7]
    alt_min = res.x[:, 2].min()
    print(f"{label:>22s}  miss {np.linalg.norm(res.r_f):8.3f} m  "
          f"vf {np.linalg.norm(res.v_f):7.4f} m/s  dv {dv:7.2f} m/s  "
          f"prop {m[0]-m[-1]:7.1f} kg  peak thr {res.thrust.max()/lander.t_max:5.2f}  "
          f"sat cyc {res.sat:4d}  min alt {round(alt_min, 1) + 0.0:7.1f} m")


def main():
    L = Lander()
    print("Initial state r0 =", R0, "m, v0 =", V0, "m/s,  m0 =", L.m0, "kg")
    print(f"Thrust/weight at high gate: {L.t_max/(L.m0*1.62):.2f} lunar g\n")

    print("Time-to-go trade (unconstrained energy-optimal arc)")
    print("   Gamma      T [s]   J_E/T [m^2/s^4]   ideal |u| peak [m/s^2]   ideal dv [m/s]")
    for gamma in [0.05, 0.2, 0.5, 1.0, 2.0, 5.0, 20.0]:
        T = optimal_tgo(R0, V0, gamma)
        a0, a1 = profile_coeffs(R0, V0, T)
        tau = np.linspace(0, T, 4001)
        u = np.linalg.norm(a0[None] + a1[None] * tau[:, None], axis=1)
        dv = np.trapezoid(u, tau)
        print(f"  {gamma:6.2f}  {T:8.2f}  {energy(R0, V0, T)/T:14.3f}  "
              f"{u.max():18.3f}  {dv:18.2f}")

    gamma = 0.5
    T = optimal_tgo(R0, V0, gamma)
    print(f"\nNominal: Gamma = {gamma}, T = {T:.2f} s")
    nominal = simulate(R0, V0, T, L)
    summary("10 Hz, limited", nominal, L)
    summary("1 Hz, limited", simulate(R0, V0, T, L, dt_guid=1.0, substeps=50), L)
    summary("10 Hz, no limits", simulate(R0, V0, T, L, limits=False), L)
    m = nominal.x[:, 6]
    print(f"  rocket-equation check: m0 exp(-dv/ve) = "
          f"{L.m0*np.exp(-nominal.x[-1,7]/(L.isp*G0)):.4f} kg vs m_f = {m[-1]:.4f} kg")

    print("\nDelta-v vs T on the ideal arc, against the full-thrust limit "
          f"{L.t_max/L.m0:.3f} m/s^2 at high-gate mass")
    for T_try in [44, 48, 52, 56, 60, 64, 70, 80, 100]:
        pk = ideal_peak_accel(R0, V0, T_try)
        print(f"  T = {T_try:5.1f} s  dv {ideal_dv(R0, V0, T_try):7.2f} m/s  "
              f"peak |u| {pk:5.3f}  {'ok' if pk <= L.t_max / L.m0 else 'SATURATES'}")
    Tf = feasible_tgo_window(R0, V0, L)
    gam_f = ((2 * V0 @ V0 * Tf**2 + 12 * R0 @ V0 * Tf + 18 * R0 @ R0) / Tf**4
             - 0.5 * G_MOON @ G_MOON)
    print(f"Shortest feasible T = {Tf:.2f} s (equivalent Gamma = {gam_f:.3f}); "
          "dv rises with T past the saturated region, so this edge is the\n"
          "lowest-dv time-to-go this guidance law can fly without saturating")
    best = simulate(R0, V0, Tf, L)
    summary("best feasible T", best, L)
    print(f"  saving vs Gamma = {gamma}: {nominal.x[-1,7]-best.x[-1,7]:.2f} m/s, "
          f"{(nominal.x[0,6]-nominal.x[-1,6])-(best.x[0,6]-best.x[-1,6]):.1f} kg")

    print("\nShortening T until the engine saturates (10 Hz, limited)")
    for T_try in [T, 0.8 * T, 0.65 * T, 0.55 * T, 0.5 * T]:
        summary(f"T = {T_try:6.2f} s", simulate(R0, V0, T_try, L), L)

    print("\nStretching T: the min-energy arc hovers and spends propellant on gravity")
    for T_try in [1.5 * T, 2.0 * T, 3.0 * T]:
        summary(f"T = {T_try:6.2f} s", simulate(R0, V0, T_try, L), L)


if __name__ == "__main__":
    main()
