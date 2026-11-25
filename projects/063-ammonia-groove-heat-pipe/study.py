# Astron Built this AP
# Project 063 | Thermal analysis | Axially grooved aluminium/ammonia heat pipe transport limits
# Author: Aditya (aditya-astron)
# Date: 2026-11-25 IST
"""
Size the grooves of a 12.7 mm aluminium/ammonia heat pipe that carries an
80 W avionics box to a radiator (Le 0.15 m, La 0.60 m, Lc 0.30 m) in orbit,
with 1.5x margin from -20 to +60 C, and that can still be acceptance tested
on the ground at 2.54 mm adverse tilt.

    python study.py
"""
import numpy as np

import heatpipe as hp

C0 = 273.15
Q_REQ = 80.0
MARGIN = 1.5
T_OPS = np.arange(-20.0, 61.0, 5.0) + C0
TILT_MM = 2.54          # standard 0.1 in adverse elevation for 1-g tests
Q_GROUND_MIN = 40.0     # must be well above rig losses to mean anything


def worst_case(p):
    """Minimum over the operating range of the lowest 0-g limit."""
    return min(min(p.limits(T).values()) for T in T_OPS)


def ground_q(p, T=20.0 + C0, tilt=TILT_MM):
    return p.q_capillary(hp.ammonia(T), tilt_mm=tilt, one_g=True)


def section(title):
    print("\n" + title)
    print("-" * len(title))


def main():
    section("1. Groove hydraulics: exact series vs finite differences")
    print("  w [mm]  d [mm]   series C [m^4]   FD (Richardson)   rel diff")
    for w, d in [(0.5, 1.2), (1.0, 1.0), (1.4, 0.4)]:
        a = hp.groove_flow_coeff(w * 1e-3, d * 1e-3)
        b = hp.groove_flow_coeff_fd(w * 1e-3, d * 1e-3, 80)
        print(f"  {w:5.2f}  {d:5.2f}   {a:.6e}     {b:.6e}     {b / a - 1:+.1e}")

    base = hp.Pipe()
    section(f"2. Baseline w = 0.50 mm, d = 1.20 mm, N = {base.N}, "
            f"Leff = {base.Leff:.3f} m")
    print("  T [C]  merit [W/m^2]  cap 0g [W]  1g+2.54mm [W]  entrain   sonic"
          "      viscous   Re_v    Ma_v")
    for Tc in range(-40, 81, 20):
        s = hp.ammonia(Tc + C0)
        L = base.limits(Tc + C0)
        Re, Ma = base.vapour_state(s, L["capillary_0g"])
        print(f"  {Tc:4d}   {s.merit:.3e}     {L['capillary_0g']:7.1f}     "
              f"{ground_q(base, Tc + C0):7.1f}      {L['entrainment']:7.0f}  "
              f"{L['sonic']:9.3g}  {L['viscous']:9.3g}  {Re:5.0f}  {Ma:.1e}")
    s20 = hp.ammonia(20 + C0)
    rl, rv = base.liquid_resistance(s20), base.vapour_resistance(s20)
    print(f"  at 20 C the vapour core takes {100 * rv / (rl + rv):.1f} % of the "
          f"viscous drop; capillary head 2 sigma/w = {base.cap_pressure(s20):.1f} Pa, "
          f"1-g circumferential head rho g Dv = {base.gravity_head(s20, 0, True):.1f} Pa")

    section("3. Effective length from the axial pressure march (20 C, 100 W)")
    Q = 100.0
    cases = {"uniform evaporator": hp.uniform_profile(base, Q)}
    for name, (a, b) in {"box on outboard 3 cm": (0.0, 0.03),
                         "box on inboard 3 cm": (base.Le - 0.03, base.Le)}.items():
        def q(x, a=a, b=b):
            if a <= x < b:
                return Q / (b - a)
            if x >= base.Le + base.La:
                return -Q / base.Lc
            return 0.0
        cases[name] = q
    R = (rl + rv)
    for name, q in cases.items():
        x, Qx, dP = hp.march(base, s20, q, Q)
        leff = dP[0] / (R * Q)
        print(f"  {name:22s} dP at end cap {dP[0]:6.2f} Pa  ->  Leff = {leff:.4f} m"
              f"  (energy closes to {Qx[-1]:+.1e} W)")
    print(f"  closed form La + (Le + Lc)/2 = {base.Leff:.4f} m; "
          f"3 cm boxes, closed form {base.La + base.Le - 0.015 + base.Lc / 2:.4f}"
          f" / {base.La + 0.015 + base.Lc / 2:.4f} m")

    section("4. Groove width trade (d = 1.6 mm, deepest the 0.5 mm wall allows)")
    print("  w [mm]   N   worst 0g [W]  binding       1g 20C 2.54mm [W]  "
          "static height 20C [mm]")
    for w in np.arange(0.30, 1.51, 0.10):
        p = hp.Pipe(w=w * 1e-3, d=1.6e-3)
        s = s20
        h_static = (p.cap_pressure(s) / (s.rho_l * hp.G0) - p.Dv) * 1e3
        Tw = min(T_OPS, key=lambda T: min(p.limits(T).values()))
        Lw = p.limits(Tw)
        bind = f"{min(Lw, key=Lw.get).split('_')[0]} {Tw - C0:+.0f}C"
        print(f"  {w:5.2f}   {p.N:2d}    {worst_case(p):7.1f}     {bind:16s}"
              f"{ground_q(p):7.1f}            {h_static:6.1f}")

    section("5. Design selection")
    feas = []
    for w in np.arange(0.25, 1.51, 0.05):
        for d in np.arange(0.4, 1.61, 0.1):
            p = hp.Pipe(w=w * 1e-3, d=d * 1e-3)
            q0, qg = worst_case(p), ground_q(p)
            if q0 >= MARGIN * Q_REQ and qg >= Q_GROUND_MIN:
                feas.append((q0, qg, p))
    print(f"  {len(feas)} of the 26 x 13 grid meet 0-g >= {MARGIN * Q_REQ:.0f} W "
          f"(-20..60 C) and 1-g >= {Q_GROUND_MIN:.0f} W at 2.54 mm")
    # Among feasible designs take the one with the best 0-g capacity; the
    # ground-test requirement is what keeps the width down.
    q0, qg, best = max(feas, key=lambda f: f[0])
    print(f"  selected w = {best.w * 1e3:.2f} mm, d = {best.d * 1e3:.1f} mm, "
          f"N = {best.N}: worst 0-g {q0:.0f} W ({q0 / Q_REQ:.2f}x), "
          f"ground {qg:.0f} W")
    worst_T = min(T_OPS, key=lambda T: min(best.limits(T).values()))
    print(f"  worst case at {worst_T - C0:.0f} C, binding limit: "
          f"{min(best.limits(worst_T), key=best.limits(worst_T).get)}")

    section("6. Ground test: Q_max vs adverse tilt at 20 C, extrapolated to 0 g")
    # Acceptance tests measure dryout power at a few small adverse tilts.
    # Q_max is linear in tilt (capillary-limited, laminar), and the 1-g
    # circumferential head is worth exactly Dv of extra tilt, so the 0-g
    # capacity is the fit evaluated at tilt = -Dv.
    tilts = np.array([0.0, 1.27, 2.54, 3.81])
    qs = np.array([ground_q(best, tilt=t) for t in tilts])
    for t, qv in zip(tilts, qs):
        print(f"  tilt {t:5.2f} mm   Q_max {qv:6.1f} W")
    X = np.column_stack([tilts, np.ones_like(tilts)])
    coef, *_ = np.linalg.lstsq(X, qs, rcond=None)
    x0 = np.array([-best.Dv * 1e3, 1.0])
    q_extrap = x0 @ coef
    print(f"  fit {coef[0]:.2f} W/mm, zero-flow tilt {-coef[1] / coef[0]:.2f} mm; "
          f"at tilt = -Dv: {q_extrap:.1f} W (model 0-g at 20 C: "
          f"{best.q_capillary(s20):.1f} W)")
    q0g = best.q_capillary(s20)
    Re_t = best.vapour_state(s20, qs[0])[0]
    Re_0 = best.vapour_state(s20, q0g)[0]
    print(f"  the {100 * (q_extrap / q0g - 1):+.1f} % bias is the vapour core: Re "
          f"{Re_t:.0f} in the tests, {Re_0:.0f} (turbulent) at the 0-g limit")
    # The extrapolation lever arm amplifies measurement scatter.
    gain = np.sqrt(x0 @ np.linalg.inv(X.T @ X) @ x0)
    print(f"  noise gain at -Dv: {gain:.2f} (vs {np.sqrt(1 / len(tilts)):.2f} "
          f"inside the data); +/-3 W per dryout point -> +/-{3 * gain:.0f} W "
          f"on the 0-g claim")

if __name__ == "__main__":
    main()
