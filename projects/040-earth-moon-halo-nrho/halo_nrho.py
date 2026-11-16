# Astron Built this AP
# Project 040 | Orbital mechanics | Earth-Moon L2 halo family to the 9:2 NRHO
# Author: Aditya (aditya-astron)
# Date: 2026-11-16 IST
"""
Walk the Earth-Moon L2 southern halo family from its bifurcation off the
planar Lyapunov family down to the 9:2 synodic-resonant near-rectilinear
halo orbit (period = synodic month / 4.5) and report its geometry and
stability.

Run: python3 halo_nrho.py
"""

import numpy as np
from scipy.optimize import brentq

import cr3bp as c
import families as f
import richardson as r

DAY = 86400.0
MOON = np.array([1 - c.MU, 0.0, 0.0])
T_NRHO_92 = c.SYNODIC_MONTH / 4.5 / c.T_STAR


def lyapunov_orbit(x0, vy_guess):
    """Planar L2 Lyapunov orbit through (x0, 0) on the Moon side."""
    s = np.array([x0, 0, 0, 0, vy_guess, 0.0])
    return f.correct(s, free=(f.IVY,), targets=(3,))


def lyapunov_family(n_steps=60, ds=0.01):
    g1, g2, x1, x2 = c.libration_points()
    rc = r.coefficients(2)
    Ax = 0.005                       # linear seed, nondimensional
    seed = np.array([x2 - Ax, 0, 0, 0, rc.k * rc.lam * Ax, 0.0])
    o = f.correct(seed, free=(f.IVY,), targets=(3,))
    return f.continue_family(o, ds, n_steps, free=(f.IX, f.IVY), targets=(3,),
                             direction=np.array([-1.0, 0.0]))


def halo_bifurcation(lyap):
    """Locate the Lyapunov orbit whose out-of-plane pair sits at +1."""
    tr = np.array([f.lyapunov_out_of_plane_trace(o) for o in lyap])
    i = int(np.argmax(tr > 2.0))
    lo, hi = lyap[i - 1], lyap[i]

    def g(x0):
        vy = np.interp(x0, [hi.s0[0], lo.s0[0]], [hi.s0[4], lo.s0[4]])
        return f.lyapunov_out_of_plane_trace(lyapunov_orbit(x0, vy)) - 2.0

    x_b = brentq(g, hi.s0[0], lo.s0[0], xtol=1e-13)
    vy = np.interp(x_b, [hi.s0[0], lo.s0[0]], [hi.s0[4], lo.s0[4]])
    return lyapunov_orbit(x_b, vy)


def halo_from_richardson(Az_km):
    """Third-order guess on the far side (tau1 = pi), corrected at fixed z0."""
    rc = r.coefficients(2)
    guess, T_guess = r.halo_state(rc, Az_km, c.L_STAR, np.pi, south=True)
    guess[[1, 3, 5]] = 0.0           # exact perpendicular crossing
    orb = f.correct(guess, free=(f.IX, f.IVY))
    return guess, T_guess, orb


def perilune_apolune(orb, n=20001):
    _, sol = c.propagate(orb.s0, orb.period, dense=True)
    X = sol.sol(np.linspace(0, orb.period, n))
    d = np.linalg.norm(X[:3].T - MOON, axis=1)
    # Both apses of these symmetric orbits sit on the x-z plane crossings.
    rp = min(np.linalg.norm(orb.s_half[:3] - MOON), d.min())
    ra = max(np.linalg.norm(orb.s0[:3] - MOON), d.max())
    return rp * c.L_STAR, ra * c.L_STAR


def walk_to_nrho(start, ds=0.01, max_steps=400):
    fam = f.continue_family(start, ds, max_steps, direction=np.array([0, -1.0, 0]),
                            stop=lambda q: q.period < T_NRHO_92)
    nrho = f.correct_fixed_period(fam[-1].s0, T_NRHO_92)
    return fam, nrho


def main():
    g1, g2, x1, x2 = c.libration_points()
    print("Earth-Moon CR3BP")
    print(f"  mu = {c.MU:.9f}, L* = {c.L_STAR:.0f} km, T* = {c.T_STAR / DAY:.5f} d")
    print(f"  L1 x = {x1:.6f} (gamma {g1:.6f}),  L2 x = {x2:.6f} (gamma {g2:.6f})")
    print(f"  C(L1) = {c.jacobi(np.array([x1, 0, 0, 0, 0, 0])):.6f},"
          f"  C(L2) = {c.jacobi(np.array([x2, 0, 0, 0, 0, 0])):.6f}")

    rc = r.coefficients(2)
    ax_min = np.sqrt(-rc.delta / rc.l1) * rc.gamma * c.L_STAR
    print("\nRichardson third-order, L2")
    print(f"  c2 {rc.c2:.6f}  c3 {rc.c3:.6f}  c4 {rc.c4:.6f}  lambda {rc.lam:.6f}  k {rc.k:.6f}")
    print(f"  l1 {rc.l1:.5f}  l2 {rc.l2:.5f}  Delta {rc.delta:.5f}  -> Ax_min {ax_min:.0f} km")

    lyap = lyapunov_family()
    bif = halo_bifurcation(lyap)
    ax_bif = (bif.s_half[0] - bif.s0[0]) / 2 * c.L_STAR
    print("\nLyapunov -> halo bifurcation (out-of-plane pair at +1)")
    print(f"  C = {bif.C:.6f}, T = {bif.period * c.T_STAR / DAY:.4f} d,"
          f" Ax = {ax_bif:.0f} km (Richardson Ax_min {ax_min:.0f} km)")

    Az0 = 10000.0
    guess, T_guess, h0 = halo_from_richardson(Az0)
    dx = (h0.s0 - guess) * c.L_STAR
    print(f"\nRichardson guess at Az = {Az0:.0f} km, corrected at fixed z0")
    print(f"  dx0 = {dx[0]:.1f} km, dvy0 = {dx[4] * 1000 / c.T_STAR:.2f} m/s,"
          f" T {T_guess * c.T_STAR / DAY:.4f} -> {h0.period * c.T_STAR / DAY:.4f} d")

    fam, nrho = walk_to_nrho(h0)
    print(f"\nPseudo-arclength continuation ({len(fam)} orbits)")
    print("   z0 [km]   rp [km]   ra [km]   T [d]      C        nu1      nu2")
    for q in fam[::5] + [nrho]:
        rp, ra = perilune_apolune(q, 4001)
        _, nus = f.stability(q)
        print(f"  {q.s0[2] * c.L_STAR:8.0f} {rp:9.0f} {ra:9.0f} {q.period * c.T_STAR / DAY:7.3f}"
              f" {q.C:9.5f} {nus[0]:8.3f} {nus[1]:7.3f}")

    rp, ra = perilune_apolune(nrho)
    ev, nus = f.stability(nrho)
    end, _ = c.propagate(nrho.s0, nrho.period)
    print("\n9:2 NRHO (southern, perilune over the north pole)")
    print(f"  s0 = [{', '.join(f'{v:.10f}' for v in nrho.s0)}]")
    print(f"  T = {nrho.period:.8f} = {nrho.period * c.T_STAR / DAY:.5f} d,  C = {nrho.C:.8f}")
    print(f"  perilune radius {rp:.0f} km (altitude {rp - c.R_MOON:.0f} km), apolune radius {ra:.0f} km")
    print(f"  stability indices {nus[0]:.4f}, {nus[1]:.4f};  max |lambda| {np.abs(ev).max():.4f}")
    print(f"  one-period closure error {np.linalg.norm(end[:3] - nrho.s0[:3]) * c.L_STAR * 1e6:.1f} mm")


if __name__ == "__main__":
    main()
