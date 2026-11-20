# Astron Built this AP
# Project 049 | Orbital mechanics | Conjunction collision probability and avoidance burn sizing
# Author: Aditya (aditya-astron)
# Date: 2026-11-20 IST
"""
A 700 km sun-synchronous satellite and a tracked debris fragment at TCA.

1. Build both states and RTN covariances, combine them in ECI and project to
   the encounter plane.
2. Compute Pc four ways and confirm it with a 3D Monte Carlo that samples
   position and velocity error.
3. Check the short-term encounter assumption: encounter duration and the
   curvature error accumulated across it.
4. Covariance scaling: where the actual covariance sits relative to the
   maximum-Pc scale (is a low Pc trustworthy or just diluted?).
5. Size the smallest along-track burn, at several lead times, that drives Pc
   below the 1e-4 manoeuvre threshold; CW prediction checked against two-body
   propagation of the burn.
"""

import math

import numpy as np
from scipy.optimize import brentq

import conjunction as cj

HBR = 15.0            # combined hard-body radius, m (10 m sat + 5 m fragment)
PC_THRESHOLD = 1e-4   # typical manoeuvre trigger
ALT = 700e3
INC = math.radians(98.19)

# 1-sigma RTN position uncertainties at TCA, m (along-track dominates)
SIG_PRIMARY = np.array([25.0, 180.0, 15.0])
SIG_SECONDARY = np.array([60.0, 420.0, 45.0])
SIG_VEL = 0.15        # m/s, per axis, both objects (used only in Monte Carlo)

CROSSING_DEG = 72.0
MISS_RTN = np.array([42.0, 0.0, -78.0])  # secondary minus primary, primary RTN, m


def scenario():
    a = cj.R_EARTH + ALT
    vc = math.sqrt(cj.MU_EARTH / a)
    raan, u = math.radians(35.0), math.radians(51.0)
    # Primary on a circular orbit, argument of latitude u.
    cr, sr = math.cos(raan), math.sin(raan)
    ci, si = math.cos(INC), math.sin(INC)
    p_hat = np.array([cr, sr, 0.0])
    q_hat = np.array([-sr * ci, cr * ci, si])
    r1 = a * (math.cos(u) * p_hat + math.sin(u) * q_hat)
    v1 = vc * (-math.sin(u) * p_hat + math.cos(u) * q_hat)
    q1 = cj.rtn_basis(r1, v1)

    # Secondary: same radius, velocity rotated about the radial by the
    # crossing angle (a fragment on a different plane through the same point).
    ang = math.radians(CROSSING_DEG)
    t_hat, n_hat = q1[:, 1], q1[:, 2]
    v2 = vc * (math.cos(ang) * t_hat + math.sin(ang) * n_hat)
    r2 = r1 + q1 @ MISS_RTN
    # Remove any offset along the relative velocity so r1, r2 are truly at TCA.
    dv = v2 - v1
    dr = r2 - r1
    r2 = r2 - (dr @ dv) / (dv @ dv) * dv
    c1 = cj.rtn_to_eci_cov(r1, v1, np.diag(SIG_PRIMARY ** 2))
    c2 = cj.rtn_to_eci_cov(r2, v2, np.diag(SIG_SECONDARY ** 2))
    return r1, v1, r2, v2, c1, c2


def main():
    rng = np.random.default_rng(20261120)
    r1, v1, r2, v2, c1, c2 = scenario()
    dr, dv = r2 - r1, v2 - v1
    cov = c1 + c2
    m, c = cj.project_encounter(dr, dv, cov)
    (sx, sy), _ = cj._principal(m, c)
    t_ca, d_ca = cj.rectilinear_miss(dr, dv)

    print("Conjunction at TCA")
    print(f"  altitude                 {ALT / 1e3:.0f} km, crossing angle {CROSSING_DEG:.0f} deg")
    print(f"  relative speed           {np.linalg.norm(dv):.1f} m/s")
    print(f"  miss distance            {np.linalg.norm(dr):.2f} m (rectilinear check {d_ca:.2f} m, t = {t_ca:.1e} s)")
    print(f"  projected sigmas         {sx:.1f} m x {sy:.1f} m (aspect {sx / sy:.2f})")
    print(f"  Mahalanobis miss         {math.sqrt(m @ np.linalg.solve(c, m)):.3f}")
    print(f"  hard-body radius         {HBR:.1f} m")

    print("\nPc by method")
    pf = cj.pc_foster(m, c, HBR)
    p1 = cj.pc_1d(m, c, HBR)
    pch = cj.pc_chan(m, c, HBR)
    psm = cj.pc_small(m, c, HBR)
    print(f"  Foster polar quadrature  {pf:.6e}")
    print(f"  Alfano 1D integral       {p1:.6e}   rel diff {abs(p1 - pf) / pf:.1e}")
    print(f"  Chan equivalent-area     {pch:.6e}   rel err {(pch - pf) / pf:+.2e}")
    print(f"  small-object             {psm:.6e}   rel err {(psm - pf) / pf:+.2e}")

    n_mc = 4_000_000
    cov_v = np.eye(3) * 2.0 * SIG_VEL ** 2
    p_mc, hits = cj.pc_monte_carlo(dr, dv, cov, cov_v, HBR, n_mc, rng)
    se = math.sqrt(pf * (1.0 - pf) / n_mc)
    print(f"  Monte Carlo 3D (n={n_mc:,})  {p_mc:.6e}   {hits} hits, "
          f"z = {(p_mc - pf) / se:+.2f}")

    print("\nShort-term encounter assumption")
    # Time for the error ellipsoid to pass: 5 sigma along the relative velocity.
    sig_z = math.sqrt(float(dv @ cov @ dv) / (dv @ dv))
    t_half = 5.0 * sig_z / np.linalg.norm(dv)
    da = -cj.MU_EARTH * (r2 / np.linalg.norm(r2) ** 3 - r1 / np.linalg.norm(r1) ** 3)
    curve = 0.5 * np.linalg.norm(da) * t_half ** 2
    print(f"  +/-5 sigma transit time  {2 * t_half * 1e3:.1f} ms")
    print(f"  relative accel           {np.linalg.norm(da):.2e} m/s^2")
    print(f"  curvature over transit   {curve:.1e} m  (vs HBR {HBR:.0f} m)")

    print("\nCovariance scaling C -> k^2 C")
    k_star, pmax = cj.max_pc_numeric(m, c, HBR)
    k_an, pmax_an = cj.max_pc_small_object(m, c, HBR)
    print(f"  max Pc numeric           {pmax:.4e} at k = {k_star:.4f}")
    print(f"  small-object closed form {pmax_an:.4e} at k = {k_an:.4f}")
    print("  k      Pc")
    for k in sorted((0.25, 0.5, k_star, 1.0, 2.0, 4.0, 8.0)):
        print(f"  {k:5.3f}  {cj.pc_scaled(m, c, HBR, k):.3e}")
    side = "below" if 1.0 < k_star else "above"
    print(f"  actual covariance (k=1) is {side} the peak scale: "
          + ("growing uncertainty raises Pc (not diluted)" if side == "below"
             else "Pc is diluted; a lower Pc here is not reassurance"))

    print(f"\nAvoidance: smallest along-track burn for Pc < {PC_THRESHOLD:.0e}")
    a = np.linalg.norm(r1)
    n = math.sqrt(cj.MU_EARTH / a ** 3)
    period = 2.0 * math.pi / n
    q1 = cj.rtn_basis(r1, v1)
    e = cj.encounter_frame(dr, dv)

    def pc_after(dvt, lead):
        # Manoeuvre moves the primary by d; secondary minus primary shifts by -d.
        d_eci = q1 @ cj.cw_impulse_displacement(n, lead, (0.0, dvt, 0.0))
        m_new = (e @ (dr - d_eci))[:2]
        return cj.pc_foster(m_new, c, HBR)

    print("  lead [orbits]  dv [mm/s]   along-track shift [m]  CW vs two-body [m]  Pc after")
    for lead_orb in (0.5, 1.0, 1.5, 2.0, 5.0):
        lead = lead_orb * period
        best = None
        for sign in (+1.0, -1.0):
            grid = sign * np.geomspace(1e-5, 1.0, 400)
            vals = np.array([pc_after(g, lead) for g in grid])
            idx = np.nonzero(vals < PC_THRESHOLD)[0]
            if len(idx) == 0:
                continue
            j = idx[0]
            lo = grid[j - 1] if j > 0 else 0.0
            root = brentq(lambda x: pc_after(x, lead) - PC_THRESHOLD, lo, grid[j],
                          xtol=1e-12)
            if best is None or abs(root) < abs(best):
                best = root
        d_cw = cj.cw_impulse_displacement(n, lead, (0.0, best, 0.0))
        d_tb = cj.twobody_impulse_displacement(r1, v1, lead, (0.0, best, 0.0))
        d_tb_eci = q1 @ d_tb
        m_tb = (e @ (dr - d_tb_eci))[:2]
        print(f"  {lead_orb:13.1f}  {best * 1e3:+9.2f}   {d_cw[1]:+20.1f}  "
              f"{np.linalg.norm(d_cw - d_tb):18.3f}  {cj.pc_foster(m_tb, c, HBR):.3e}")


if __name__ == "__main__":
    main()
