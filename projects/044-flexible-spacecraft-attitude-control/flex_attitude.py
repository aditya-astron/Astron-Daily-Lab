# Astron Built this AP
# Project 044 | PID and control systems | Attitude control of a spacecraft with flexible solar arrays
# Author: Aditya (aditya-astron)
# Date: 2026-11-17 IST
"""
A 2 Hz PD attitude loop on a spacecraft whose two solar panels have their
first antisymmetric mode at 0.4 Hz. Compare five ways of living with the
mode: do nothing, notch it, notch below it, roll the loop off, or sample
faster. Each is judged nominally and over a +/- modal-frequency error, since
on-orbit panel frequencies are rarely known to better than 10-20 %.

Run: python3 flex_attitude.py
"""

import numpy as np
from scipy.optimize import brentq

import control as c
from flexsat import FlexSat

WN, ZETA = 0.2, 0.7          # rigid-body target: wn rad/s, damping
N_MODES = 3                  # modal effective inertia > 99.9 %
STEP = np.radians(0.5)       # small slew used for the time response


def designs(sat):
    wp = sat.free_free_fe(1)[0]
    return [
        ("PD, 2 Hz", 0.5, None),
        ("PD + notch at w_p, 2 Hz", 0.5, c.notch(wp, 0.5, 0.7)),
        ("PD + notch at 0.8 w_p, 2 Hz", 0.5, c.notch(0.8 * wp, 0.5, 0.7)),
        ("PD + low-pass 0.8 rad/s, 2 Hz", 0.5, c.lowpass2(0.8, 0.5, 0.5)),
        ("PD, 10 Hz", 0.1, None),
    ]


def loop(sat, T, filt, **kw):
    return c.Loop(sat, T, WN, ZETA, filt, n_modes=kw.pop("n_modes", N_MODES), **kw)


def least_damped_flex(lp):
    """Least-damped closed-loop pole above 1 rad/s (panel modes, incl. aliased ones)."""
    return min((z, w) for w, z in lp.mode_table() if w > 1.0)


def edge(f, x0, step, limit):
    """First root of f moving from x0 (f(x0) < 0) in steps of `step`, or limit."""
    x = x0
    while (x + step - limit) * np.sign(step) <= 0:
        if f(x + step) >= 0:
            return brentq(f, x, x + step, xtol=1e-6)
        x += step
    return limit


def freq_band(sat, T, filt, lo=0.3, hi=2.5):
    f = lambda s: loop(sat, T, filt, freq_scale=s).spectral_radius() - 1.0
    if f(1.0) >= 0:
        return None
    return edge(f, 1.0, -0.01, lo), edge(f, 1.0, 0.01, hi)


def gain_band(sat, T, filt, hi=100.0):
    f = lambda g: loop(sat, T, filt, gain=g).spectral_radius() - 1.0
    if f(1.0) >= 0:
        return None
    lo_edge = edge(lambda x: f(np.exp(x)), 0.0, -0.05, np.log(1e-3))
    hi_edge = edge(lambda x: f(np.exp(x)), 0.0, 0.05, np.log(hi))
    return np.exp(lo_edge), np.exp(hi_edge)


def slew(sat, T, filt, t_end=150.0):
    lp = loop(sat, T, filt)
    N = int(round(t_end / T))
    X = lp.simulate(STEP, N)
    t = np.arange(N + 1) * T
    th = X[:, 0]
    u = X[:, -1]
    tip = sat.tip_deflection(X[:, :2 * N_MODES + 2], N_MODES)
    out = abs(th - STEP) > 0.02 * STEP
    ts = t[np.nonzero(out)[0][-1] + 1] if out.any() and not out[-1] else np.inf
    return dict(t=t, theta=th, u=u, tip=tip, ts=ts,
                u_pk=np.max(np.abs(u)), tip_pk=np.max(np.abs(tip)),
                tip_end=np.max(np.abs(tip[t > t_end - 20])))


def critical_period(sat, scale, n_modes=1):
    f = lambda T: loop(sat, T, None, freq_scale=scale, n_modes=n_modes).spectral_radius() - 1
    Ts = np.linspace(0.02, 1.0, 50)
    v = [f(T) for T in Ts]
    i = next(i for i in range(len(v)) if v[i] > 0)
    return brentq(f, Ts[i - 1], Ts[i], xtol=1e-7)


def phase_rule_period(wp, wn=WN, zeta=ZETA):
    """Rate feedback damps a collocated mode while the PD lead beats the
    ZOH + one-sample lag: 1.5 w T < atan(2 zeta w / wn)."""
    return np.arctan(2 * zeta * wp / wn) / (1.5 * wp)


def main():
    sat = FlexSat()
    wp = sat.free_free_fe(3)
    print("Structure")
    print(f"  J_t = {sat.J_t:.0f} kg m^2 (hub {sat.J_hub:.0f}, panels {2 * sat.J_panel:.0f})")
    for k in range(3):
        print(f"  mode {k + 1}: cantilever {sat.w_c[k]:7.3f} rad/s, L = {sat.L[k]:6.2f} kg^0.5 m,"
              f" free-free {wp[k]:7.3f} rad/s ({wp[k] / 2 / np.pi:.3f} Hz)")
    eff = np.cumsum(sat.L**2) / (2 * sat.J_panel)
    print(f"  modal effective inertia (panels): 1 mode {eff[0]:.4f}, 2 {eff[1]:.4f}, 3 {eff[2]:.4f}")
    print(f"  one-mode estimate w_c / sqrt(1 - L1^2/J_t) = "
          f"{sat.w_c[0] / np.sqrt(1 - sat.L[0]**2 / sat.J_t):.4f} rad/s vs FE {wp[0]:.4f}")
    print(f"  at 2 Hz, modes 2 and 3 fold to {abs(wp[1] - 4 * np.pi):.2f} and "
          f"{abs(wp[2] - 8 * np.pi):.2f} rad/s\n")

    print(f"Rigid-body PD: wn = {WN} rad/s, zeta = {ZETA}, Kp = {sat.J_t * WN**2:.1f} N m/rad,"
          f" Kd = {2 * ZETA * WN * sat.J_t:.1f} N m s/rad\n")

    hdr = (f"{'design':32s} {'stable':>6s} {'worst flex pole':>17s} {'rigid PM':>8s}"
           f" {'GM':>7s} {'w_p range (x nominal)':>22s}")
    print(hdr)
    print("-" * len(hdr))
    rows = []
    for name, T, filt in designs(sat):
        lp = loop(sat, T, filt)
        stable = lp.spectral_radius() < 1
        zf, wf = least_damped_flex(lp)
        _, gc = lp.margins()
        pm = min(p for w, p in gc if w < 1.0)
        gb = gain_band(sat, T, filt)
        fb = freq_band(sat, T, filt)
        gtxt = f"{20 * np.log10(gb[1]):.1f}dB" if gb else "-"
        ftxt = f"{fb[0]:.2f} - {fb[1]:.2f}" if fb else "-"
        print(f"{name:32s} {'yes' if stable else 'NO':>6s} {zf:+7.4f} @{wf:5.2f} r/s"
              f" {pm:7.1f}d {gtxt:>7s} {ftxt:>22s}")
        rows.append((name, T, filt, stable))
    print("  worst flex pole: damping ratio and frequency of the least-damped pole above 1 rad/s"
          " (open loop 0.005)")
    print("  GM: loop-gain increase to the first unit-circle crossing (gain reductions"
          " stay stable down to 0.001x)\n")

    print(f"Slew response, {np.degrees(STEP):.1f} deg step")
    for name, T, filt, stable in rows:
        if not stable:
            r = slew(sat, T, filt, 60.0)
            print(f"  {name:32s} diverges: tip {1e3 * r['tip_pk']:.1f} mm after 60 s")
            continue
        r = slew(sat, T, filt)
        print(f"  {name:32s} settle(2%) {r['ts']:6.1f} s, peak torque {r['u_pk']:.3f} N m,"
              f" peak tip {1e3 * r['tip_pk']:.2f} mm, last 20 s {1e3 * r['tip_end']:.3f} mm")

    print("\nLongest stable sample period for plain PD vs the 1.5 w T phase rule")
    for scale in (0.6, 0.8, 1.0, 1.2, 1.5):
        Tc = critical_period(sat, scale)
        Te = phase_rule_period(scale * wp[0])
        print(f"  w_p x {scale:.1f}: exact {Tc:.4f} s, phase rule {Te:.4f} s ({100 * (Te / Tc - 1):+.1f} %)")


if __name__ == "__main__":
    main()
