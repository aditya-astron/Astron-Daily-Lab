# Astron Built this AP
# Project 003 | Orbital mechanics | Lambert targeting for a LEO rendezvous transfer
# Author: Aditya (aditya-astron)
# Date: 2026-11-01 IST
"""Two-impulse rendezvous from a 300 km parking orbit to a target in an
800 km circular orbit (coplanar). The target starts 40 deg ahead. Sweep
coast time before the burn and time of flight, solve Lambert for each pair,
and find the cheapest intercept. The answer is checked against Hohmann.
"""

import numpy as np
from scipy.optimize import minimize

from lambert import MU_EARTH, R_EARTH, lambert, propagate

ALT_CHASER = 300.0          # km
ALT_TARGET = 800.0          # km
PHASE0_DEG = 40.0           # target lead angle at t = 0


def circ_state(radius, angle, mu=MU_EARTH):
    v = np.sqrt(mu / radius)
    c, s = np.cos(angle), np.sin(angle)
    return np.array([radius * c, radius * s, 0.0]), np.array([-v * s, v * c, 0.0])


class Scenario:
    def __init__(self, alt_c=ALT_CHASER, alt_t=ALT_TARGET, phase0_deg=PHASE0_DEG, mu=MU_EARTH):
        self.mu = mu
        self.rc = R_EARTH + alt_c
        self.rt = R_EARTH + alt_t
        self.nc = np.sqrt(mu / self.rc**3)
        self.nt = np.sqrt(mu / self.rt**3)
        self.phase0 = np.radians(phase0_deg)

    def chaser(self, t):
        return circ_state(self.rc, self.nc * t, self.mu)

    def target(self, t):
        return circ_state(self.rt, self.phase0 + self.nt * t, self.mu)

    def transfer(self, t_wait, tof):
        """Total delta-v (km/s) and burns for depart at t_wait, arrive at t_wait+tof."""
        r1, vc = self.chaser(t_wait)
        r2, vt = self.target(t_wait + tof)
        v1, v2 = lambert(r1, r2, tof, self.mu)
        dv1 = np.linalg.norm(v1 - vc)
        dv2 = np.linalg.norm(vt - v2)
        return dv1 + dv2, dv1, dv2, (r1, v1, r2)

    # Hohmann reference -----------------------------------------------------
    def hohmann(self):
        """Classic Hohmann delta-v, transfer time and the coast needed to
        reach the right phase angle. Target must lead by pi - nt*tH at burn."""
        a = 0.5 * (self.rc + self.rt)
        dv1 = np.sqrt(self.mu / self.rc) * (np.sqrt(2 * self.rt / (self.rc + self.rt)) - 1)
        dv2 = np.sqrt(self.mu / self.rt) * (1 - np.sqrt(2 * self.rc / (self.rc + self.rt)))
        t_h = np.pi * np.sqrt(a**3 / self.mu)
        lead_req = np.pi - self.nt * t_h
        # Phase angle falls at (nc - nt) since the lower chaser is faster.
        t_wait = ((self.phase0 - lead_req) % (2 * np.pi)) / (self.nc - self.nt)
        return dv1 + dv2, dv1, dv2, t_h, t_wait

    def synodic(self):
        return 2 * np.pi / (self.nc - self.nt)


def safe_dv(sc, t_wait, tof):
    try:
        return sc.transfer(t_wait, tof)[0]
    except (ValueError, RuntimeError):
        return np.inf


def search(sc, n_wait=181, n_tof=121):
    """Coarse grid over one synodic period, then a local polish.

    The optimum sits next to the 180 deg singularity, so the polish uses
    Nelder-Mead (derivative-free) and the solver rejects exact 180 deg."""
    _, _, _, t_h, _ = sc.hohmann()
    waits = np.linspace(0.0, sc.synodic(), n_wait)
    tofs = np.linspace(0.3 * t_h, 1.6 * t_h, n_tof)
    grid = np.array([[safe_dv(sc, w, t) for t in tofs] for w in waits])
    i, j = np.unravel_index(np.argmin(grid), grid.shape)

    res = minimize(lambda x: safe_dv(sc, x[0], x[1]), [waits[i], tofs[j]],
                   method="Nelder-Mead",
                   options={"xatol": 1e-3, "fatol": 1e-10, "maxiter": 4000})
    return res.x[0], res.x[1], res.fun, waits, tofs, grid


def main():
    sc = Scenario()
    dvh, dvh1, dvh2, t_h, tw_h = sc.hohmann()

    print(f"Chaser   : {ALT_CHASER:.0f} km circular   period {2*np.pi/sc.nc/60:7.2f} min")
    print(f"Target   : {ALT_TARGET:.0f} km circular   period {2*np.pi/sc.nt/60:7.2f} min")
    print(f"Lead     : {PHASE0_DEG:.1f} deg at t=0     synodic {sc.synodic()/3600:6.2f} h")
    print()
    print("Hohmann reference")
    print(f"  dv1 {dvh1*1e3:7.2f} m/s  dv2 {dvh2*1e3:7.2f} m/s  total {dvh*1e3:7.2f} m/s")
    print(f"  transfer {t_h/60:6.2f} min, coast before burn {tw_h/60:7.2f} min")
    print()

    tw, tof, dv, waits, tofs, grid = search(sc)
    total, dv1, dv2, (r1, v1, r2) = sc.transfer(tw, tof)
    r_end, _ = propagate(r1, v1, tof)
    miss = np.linalg.norm(r_end - r2)

    print("Lambert optimum (grid + Nelder-Mead)")
    print(f"  coast {tw/60:7.2f} min  tof {tof/60:6.2f} min")
    print(f"  dv1 {dv1*1e3:7.2f} m/s  dv2 {dv2*1e3:7.2f} m/s  total {total*1e3:7.2f} m/s")
    print(f"  excess over Hohmann {(total-dvh)*1e3:.3f} m/s")
    print(f"  propagated miss distance {miss*1e3:.3e} m")
    print()

    # Cost of being in a hurry: best dv when arrival is forced early.
    print("Fast-intercept penalty (best coast for each tof)")
    print("   tof [min]   best dv [m/s]")
    for frac in (0.4, 0.6, 0.8, 1.0, 1.2):
        k = np.argmin(np.abs(tofs - frac * t_h))
        print(f"   {tofs[k]/60:8.2f}    {np.min(grid[:, k])*1e3:10.1f}")


if __name__ == "__main__":
    main()
