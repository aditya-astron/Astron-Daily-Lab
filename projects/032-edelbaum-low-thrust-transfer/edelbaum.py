# Astron Built this AP
# Project 032 | Orbital mechanics | Edelbaum low-thrust LEO-to-GEO transfer with plane change
# Author: Aditya (aditya-astron)
# Date: 2026-11-11 IST
"""
Low-thrust circle-to-circle transfer with a simultaneous plane change.

Edelbaum's quasi-circular model treats the orbit as a slowly shrinking or
growing circle of speed v and inclination i. Thrust acceleration f is split
by a yaw angle beta between the velocity direction and the orbit normal, and
the normal component flips sign at the antinodes (sign of cos u) so that it
always works against the inclination. Orbit-averaged:

    dv/dt = -f cos(beta)
    di/dt = -(2/pi) f sin(beta) / v

Minimising the total delta-V gives a yaw law that is constant in the
"v sin(beta)" sense, v sin(beta) = v0 sin(beta0), which integrates to

    v(V)   = sqrt(v0^2 - 2 v0 V cos(beta0) + V^2)
    beta(V) = atan2(v0 sin(beta0), v0 cos(beta0) - V)
    di(V)  = (2/pi) (beta(V) - beta0)

with V the accumulated delta-V. Solving for the end state gives

    tan(beta0) = sin(pi di/2) / (v0/v1 - cos(pi di/2))
    dV         = sqrt(v0^2 - 2 v0 v1 cos(pi di/2) + v1^2)

Because the law is written in terms of accumulated delta-V, it stays optimal
for constant thrust with a falling mass; only the clock changes.

This module provides the closed form, an independent direct-optimisation
check (piecewise-constant yaw solved by SLSQP), a full Cartesian propagation
of the steered spacecraft, and impulsive/sequential baselines.
"""

import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize, minimize_scalar

MU = 398600.4418          # km^3/s^2, Earth
G0 = 9.80665e-3           # km/s^2
R_EARTH = 6378.137        # km
R_GEO = 42164.17          # km


def circ_speed(r):
    return math.sqrt(MU / r)


# ---------------------------------------------------------------------------
# Edelbaum closed form
# ---------------------------------------------------------------------------

def edelbaum(v0, v1, di):
    """Return (dV, beta0, beta1) for a circle-to-circle transfer.

    v0, v1 in km/s, di (inclination change, magnitude) in rad.
    """
    c = math.cos(0.5 * math.pi * di)
    s = math.sin(0.5 * math.pi * di)
    dv = math.sqrt(v0 * v0 - 2.0 * v0 * v1 * c + v1 * v1)
    beta0 = math.atan2(s, v0 / v1 - c)
    beta1 = beta_of_dv(v0, beta0, dv)
    return dv, beta0, beta1


def beta_of_dv(v0, beta0, V):
    return math.atan2(v0 * math.sin(beta0), v0 * math.cos(beta0) - V)


def v_of_dv(v0, beta0, V):
    return math.sqrt(v0 * v0 - 2.0 * v0 * V * math.cos(beta0) + V * V)


def di_of_dv(v0, beta0, V):
    return 2.0 / math.pi * (beta_of_dv(v0, beta0, V) - beta0)


def averaged_check(v0, beta0, dv_total):
    """Integrate the orbit-averaged equations under the closed-form yaw law.

    Independent of the closed-form v(V) and di(V); used to confirm them.
    Independent variable is accumulated delta-V, so f drops out.
    """
    def rhs(V, y):
        v, _ = y
        b = beta_of_dv(v0, beta0, V)
        return [-math.cos(b), 2.0 / math.pi * math.sin(b) / v]

    sol = solve_ivp(rhs, (0.0, dv_total), [v0, 0.0], method="DOP853",
                    rtol=1e-12, atol=1e-14)
    return sol.y[0, -1], sol.y[1, -1]


# ---------------------------------------------------------------------------
# Direct optimisation check: piecewise-constant yaw
# ---------------------------------------------------------------------------

def _logmean(a, b):
    if abs(a - b) < 1e-12 * a:
        return 0.5 * (a + b)
    return (a - b) / math.log(a / b)


def direct_optimum(v0, v1, di, n_seg=40):
    """Minimise total delta-V over n_seg constant-yaw arcs.

    On an arc of constant beta the speed is linear in V, so the inclination
    change is exact: (2/pi) sin(beta) dV / logmean(v_start, v_end).
    Decision vector: [beta_1..beta_n, dV_1..dV_n].
    """
    dv_e, b0, b1 = edelbaum(v0, v1, di)
    # Start away from the analytic answer: linear yaw ramp, equal arcs.
    x0 = np.concatenate([np.linspace(0.3, 1.2, n_seg),
                         np.full(n_seg, 1.15 * dv_e / n_seg)])

    def propagate(x):
        beta, dvs = x[:n_seg], x[n_seg:]
        v, inc = v0, 0.0
        for b, d in zip(beta, dvs):
            vn = v - math.cos(b) * d
            inc += 2.0 / math.pi * math.sin(b) * d / _logmean(v, vn)
            v = vn
        return v, inc

    cons = [{"type": "eq", "fun": lambda x: propagate(x)[0] - v1},
            {"type": "eq", "fun": lambda x: propagate(x)[1] - di}]
    bounds = [(0.0, math.pi)] * n_seg + [(1e-6, None)] * n_seg
    res = minimize(lambda x: np.sum(x[n_seg:]), x0, method="SLSQP",
                   bounds=bounds, constraints=cons,
                   options={"maxiter": 2000, "ftol": 1e-13})
    return res.fun, res.x, res.success


# ---------------------------------------------------------------------------
# Baselines
# ---------------------------------------------------------------------------

def sequential_lowthrust(v0, v1, di):
    """Spiral in-plane first, then a pure plane change at the final orbit."""
    return abs(v0 - v1) + 0.5 * math.pi * v1 * di


def hohmann_plane_change(r0, r1, di):
    """Impulsive Hohmann with the plane change optimally split across burns.

    Returns (dV_total, alpha_at_perigee).
    """
    v0, v1 = circ_speed(r0), circ_speed(r1)
    a_t = 0.5 * (r0 + r1)
    vp = math.sqrt(MU * (2.0 / r0 - 1.0 / a_t))
    va = math.sqrt(MU * (2.0 / r1 - 1.0 / a_t))

    def total(alpha):
        b1 = math.sqrt(v0 * v0 + vp * vp - 2 * v0 * vp * math.cos(alpha))
        b2 = math.sqrt(va * va + v1 * v1 - 2 * va * v1 * math.cos(di - alpha))
        return b1 + b2

    res = minimize_scalar(total, bounds=(0.0, di), method="bounded",
                          options={"xatol": 1e-12})
    return res.fun, res.x


# ---------------------------------------------------------------------------
# Cartesian propagation of the steered spacecraft
# ---------------------------------------------------------------------------

def elements(r, v):
    """Return (a, e, i) from a Cartesian state."""
    r = np.asarray(r)
    v = np.asarray(v)
    rn = np.linalg.norm(r)
    h = np.cross(r, v)
    a = 1.0 / (2.0 / rn - v @ v / MU)
    e_vec = np.cross(v, h) / MU - r / rn
    inc = math.acos(max(-1.0, min(1.0, h[2] / np.linalg.norm(h))))
    return a, np.linalg.norm(e_vec), inc


def propagate_cartesian(r0, r1, i0, i1, thrust_n, m0_kg, isp_s, rtol=1e-9):
    """Fly the Edelbaum yaw law in two-body Cartesian dynamics.

    Constant thrust, mass decreasing at mdot = T/(g0 Isp). The yaw is
    scheduled on accumulated delta-V (a state), not time, and the normal
    component switches sign with cos(u) so it always lowers |i - i1|.
    Integration stops when accumulated delta-V reaches the Edelbaum value.

    The residual eccentricity is the forced response to the in-plane thrust,
    of order f / (n^2 r); it is the size of the quasi-circular error.

    Returns a dict with the trajectory summary.
    """
    v0 = circ_speed(r0)
    v1 = circ_speed(r1)
    di = abs(i0 - i1)
    dv_tot, beta0, _ = edelbaum(v0, v1, di)
    sgn = 1.0 if i1 < i0 else -1.0

    T = thrust_n * 1e-3          # kg km/s^2
    mdot = thrust_n / (G0 * 1e3 * isp_s)

    def rhs(t, y):
        x, yy, z, vx, vy, vz, m, V = y
        rr = math.sqrt(x * x + yy * yy + z * z)
        g = -MU / rr ** 3
        # Orbit frame: velocity direction and angular-momentum direction.
        hx, hy, hz = yy * vz - z * vy, z * vx - x * vz, x * vy - yy * vx
        hn = math.sqrt(hx * hx + hy * hy + hz * hz)
        vn = math.sqrt(vx * vx + vy * vy + vz * vz)
        # Node line n = z_hat x h = (-hy, hx, 0); cos u = n.r / (|n| |r|).
        nn = math.hypot(hx, hy)
        cos_u = (-hy * x + hx * yy) / (nn * rr) if nn > 0 else 0.0
        b = math.atan2(v0 * math.sin(beta0), v0 * math.cos(beta0) - V)
        f = T / m
        ft = f * math.cos(b) / vn
        fh = -sgn * f * math.sin(b) * math.copysign(1.0, cos_u) / hn
        return [vx, vy, vz,
                g * x + ft * vx + fh * hx,
                g * yy + ft * vy + fh * hy,
                g * z + ft * vz + fh * hz,
                -mdot, f]

    def done(t, y):
        return y[7] - dv_tot
    done.terminal = True

    # Start on the ascending node of a circular orbit.
    y0 = [r0, 0.0, 0.0, 0.0, v0 * math.cos(i0), v0 * math.sin(i0),
          m0_kg, 0.0]
    t_est = m0_kg / mdot * (1.0 - math.exp(-dv_tot / (G0 * isp_s)))
    sol = solve_ivp(rhs, (0.0, 1.5 * t_est), y0, method="DOP853",
                    rtol=rtol, atol=1e-9, events=done, dense_output=False)
    yf = sol.y[:, -1]
    a, e, inc = elements(yf[:3], yf[3:6])

    return {
        "t_end": sol.t[-1],
        "t_rocket": t_est,
        "dv_target": dv_tot,
        "dv_flown": yf[7],
        "m_end": yf[6],
        "a": a, "e": e, "i": inc,
        "nfev": sol.nfev,
        "status": sol.status,
    }


# ---------------------------------------------------------------------------
# Demo
# ---------------------------------------------------------------------------

def main():
    r0 = R_EARTH + 200.0
    i0 = math.radians(28.5)
    v0, v1 = circ_speed(r0), circ_speed(R_GEO)
    di = i0

    dv, b0, b1 = edelbaum(v0, v1, di)
    print("Edelbaum LEO 200 km / 28.5 deg -> GEO / 0 deg")
    print(f"  v0 = {v0:.4f} km/s, v1 = {v1:.4f} km/s")
    print(f"  dV = {dv:.4f} km/s, beta0 = {math.degrees(b0):.3f} deg, "
          f"beta_end = {math.degrees(b1):.3f} deg")

    print("\n  Yaw and inclination schedule (accumulated dV):")
    print("     V [km/s]  r [km]    beta [deg]  i [deg]")
    for frac in np.linspace(0, 1, 6):
        V = frac * dv
        v = v_of_dv(v0, b0, V)
        print(f"    {V:8.3f}  {MU / v**2:8.0f}  {math.degrees(beta_of_dv(v0, b0, V)):9.3f}"
              f"  {math.degrees(i0 - di_of_dv(v0, b0, V)):8.3f}")

    v_end, di_end = averaged_check(v0, b0, dv)
    print(f"\n  Averaged ODE under the yaw law: v_end = {v_end:.9f} km/s, "
          f"di = {math.degrees(di_end):.9f} deg")

    dv_dir, _, ok = direct_optimum(v0, v1, di)
    print(f"  Direct optimisation (40 constant-yaw arcs, SLSQP): "
          f"dV = {dv_dir:.5f} km/s ({'converged' if ok else 'NOT converged'}), "
          f"+{(dv_dir / dv - 1) * 100:.3f} % over Edelbaum")

    seq = sequential_lowthrust(v0, v1, di)
    hoh, alpha = hohmann_plane_change(r0, R_GEO, di)
    print("\n  Baselines:")
    print(f"    Sequential low-thrust (spiral, then plane change at GEO): {seq:.4f} km/s")
    print(f"    Impulsive Hohmann, optimal split ({math.degrees(alpha):.2f} deg at perigee): "
          f"{hoh:.4f} km/s")
    print(f"    Low-thrust gravity-loss penalty vs. impulsive: {dv - hoh:.4f} km/s")

    # Electric orbit raising, constant thrust, full Cartesian propagation.
    m0, thrust, isp = 2000.0, 2.0, 1800.0
    print(f"\n  Cartesian propagation: m0 = {m0:.0f} kg, T = {thrust} N, Isp = {isp:.0f} s")
    res = propagate_cartesian(r0, R_GEO, i0, 0.0, thrust, m0, isp)
    days = res["t_end"] / 86400.0
    print(f"    flight time {days:.2f} d (rocket-equation estimate "
          f"{res['t_rocket'] / 86400:.2f} d), {res['nfev']} RHS evaluations")
    print(f"    propellant {m0 - res['m_end']:.1f} kg "
          f"(rocket equation {m0 * (1 - math.exp(-dv / (G0 * isp))):.1f} kg)")
    print(f"    final a = {res['a']:.1f} km (target {R_GEO:.1f}), "
          f"e = {res['e']:.2e}, i = {math.degrees(res['i']):.4f} deg")


if __name__ == "__main__":
    main()
