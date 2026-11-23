# Astron Built this AP
# Project 060 | Navigation / Guidance | Lunar terminal descent with energy-optimal ZEM/ZEV guidance
# Author: Aditya (aditya-astron)
# Date: 2026-11-23 IST
"""Energy-optimal (ZEM/ZEV) powered-descent guidance for a lunar lander.

Problem: bring a lander from a high-gate state a couple of kilometres up and
downrange to rest on a target pad, with a throttleable engine of limited
thrust, and pick the flight time sensibly.

Model (target-centred, flat Moon, uniform gravity; fine for the last few km):
    r'' = u + g,      g = (0, 0, -1.62) m/s^2,   u = F/m (thrust acceleration)
    m'  = -|F| / (Isp g0)

Energy-optimal control for fixed time-to-go T, minimising int |u|^2 dt with
r(T) = 0, v(T) = 0, is linear in time. In terms of the zero-effort miss and
zero-effort velocity (where you would end up if the engine were off):
    ZEM = r_f - (r + v T + g T^2 / 2)
    ZEV = v_f - (v + g T)
    u   = 6 ZEM / T^2 - 2 ZEV / T
and the energy it spends is
    J_E(T) = int |u|^2 = 12|ZEM|^2/T^3 - 12 ZEM.ZEV/T^2 + 4|ZEV|^2/T.
Re-evaluating u every cycle in closed loop reproduces the same straight line
in time (the law is its own receding-horizon solution), so any re-plan error
comes only from saturation, sampling and the stop-update logic.

Time-to-go: J_E falls monotonically with T, so free final time needs a time
penalty. Minimising  Gamma T + J_E/2  gives the stationarity condition
    (Gamma + |g|^2/2) T^4 - 2|v|^2 T^2 - 12 (r.v) T - 18 |r|^2 = 0
(r, v relative to the target, v_f = 0). The coefficient signs change once,
so by Descartes' rule there is exactly one positive root.
"""
from dataclasses import dataclass
import numpy as np

G0 = 9.80665
G_MOON = np.array([0.0, 0.0, -1.62])


def zem_zev(r, v, T, g=G_MOON, rf=np.zeros(3), vf=np.zeros(3)):
    zem = rf - (r + v * T + 0.5 * g * T * T)
    zev = vf - (v + g * T)
    return zem, zev


def accel_command(r, v, T, g=G_MOON):
    """Thrust acceleration u (excludes gravity) for time-to-go T."""
    zem, zev = zem_zev(r, v, T, g)
    return 6.0 * zem / T**2 - 2.0 * zev / T


def profile_coeffs(r, v, T, g=G_MOON):
    """u(tau) = a0 + a1 tau over the remaining arc (tau measured from now)."""
    zem, zev = zem_zev(r, v, T, g)
    a0 = 6.0 * zem / T**2 - 2.0 * zev / T
    a1 = -12.0 * zem / T**3 + 6.0 * zev / T**2
    return a0, a1


def energy(r, v, T, g=G_MOON):
    zem, zev = zem_zev(r, v, T, g)
    return 12 * zem @ zem / T**3 - 12 * zem @ zev / T**2 + 4 * zev @ zev / T


def optimal_tgo(r, v, gamma, g=G_MOON):
    """Unique positive root of the free-final-time quartic."""
    c = [gamma + 0.5 * g @ g, 0.0, -2.0 * v @ v, -12.0 * r @ v, -18.0 * r @ r]
    roots = np.roots(c)
    pos = [z.real for z in roots if abs(z.imag) < 1e-9 * abs(z) and z.real > 0]
    assert len(pos) == 1, roots
    return pos[0]


@dataclass
class Lander:
    m0: float = 12000.0        # kg at high gate
    t_max: float = 45040.0     # N, Apollo-class DPS full thrust
    throttle_min: float = 0.10 # deep throttle floor
    isp: float = 311.0         # s


@dataclass
class Result:
    t: np.ndarray
    x: np.ndarray      # rows: r(3), v(3), m, dv
    thrust: np.ndarray # commanded |F| per guidance cycle (N)
    sat: int           # guidance cycles where a thrust limit was active

    @property
    def r_f(self):
        return self.x[-1, 0:3]

    @property
    def v_f(self):
        return self.x[-1, 3:6]


def simulate(r0, v0, T0, lander=Lander(), dt_guid=0.1, substeps=10,
             t_freeze=1.0, limits=True, g=G_MOON):
    """Fly the descent with ZOH guidance at 1/dt_guid Hz and RK4 dynamics.

    Fixed final time: T counts down from T0. Below t_freeze s the command is
    no longer recomputed (u ~ 1/T^2 amplifies any residual); the last
    profile a0 + a1 tau is flown open loop, which is exact if nothing
    saturates. Mass and accumulated delta-v are integrated as states.
    """
    ve = lander.isp * G0
    fmin = lander.throttle_min * lander.t_max
    # Shrink the cycle slightly so an integer number of cycles ends exactly
    # at T0; otherwise the final-time mismatch shows up as a velocity error.
    n = int(np.ceil(T0 / dt_guid - 1e-9))
    dt_guid = T0 / n
    h = dt_guid / substeps
    x = np.concatenate([r0, v0, [lander.m0, 0.0]]).astype(float)
    ts, xs, fs = [0.0], [x.copy()], []
    sat = 0
    frozen = None
    for k in range(n):
        t = k * dt_guid
        T = T0 - t
        if T > t_freeze or frozen is None:
            a0, a1 = profile_coeffs(x[0:3], x[3:6], T, g)
            frozen = (t, a0, a1) if T <= t_freeze else None
        else:
            a0, a1 = frozen[1], frozen[2]
        t_ref = frozen[0] if frozen else t
        # Within a guidance cycle hold the command constant (flight software
        # ZOH) unless frozen, where the analytic profile is flown.
        u_hold = a0 + a1 * (t - t_ref)

        def thrust_vec(m, tau):
            u = a0 + a1 * (tau - t_ref) if frozen else u_hold
            F = m * u
            mag = np.linalg.norm(F)
            if limits:
                if mag > lander.t_max:
                    return F * lander.t_max / mag, True
                if mag < fmin:
                    d = F / mag if mag > 0 else np.array([0, 0, 1.0])
                    return d * fmin, True
            return F, False

        _, s = thrust_vec(x[6], t)
        sat += s
        fs.append(np.linalg.norm(thrust_vec(x[6], t)[0]))

        def f(tau, y):
            F, _ = thrust_vec(y[6], tau)
            mag = np.linalg.norm(F)
            return np.concatenate([y[3:6], F / y[6] + g, [-mag / ve, mag / y[6]]])

        for j in range(substeps):
            tau = t + j * h
            k1 = f(tau, x)
            k2 = f(tau + h / 2, x + h / 2 * k1)
            k3 = f(tau + h / 2, x + h / 2 * k2)
            k4 = f(tau + h, x + h * k3)
            x = x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        ts.append(t + dt_guid)
        xs.append(x.copy())
    return Result(np.array(ts), np.array(xs), np.array(fs), sat)


def ideal_peak_accel(r, v, T, g=G_MOON):
    """Peak |u| on the unconstrained arc. |a0 + a1 tau| is convex in tau,
    so the maximum over [0, T] sits at one of the two ends."""
    a0, a1 = profile_coeffs(r, v, T, g)
    return max(np.linalg.norm(a0), np.linalg.norm(a0 + a1 * T))


def ideal_dv(r, v, T, g=G_MOON, n=4001):
    a0, a1 = profile_coeffs(r, v, T, g)
    tau = np.linspace(0.0, T, n)
    return np.trapezoid(np.linalg.norm(a0[None] + a1[None] * tau[:, None], axis=1), tau)


def feasible_tgo_window(r, v, lander, g=G_MOON, t_lo=5.0, t_hi=400.0, step=0.05):
    """Shortest T whose unconstrained arc never needs more than full thrust.

    Conservative: uses the high-gate mass, so the late (lighter) part of the
    arc has extra margin. Scans, then refines the crossing by bisection.
    """
    a_max = lander.t_max / lander.m0
    Ts = np.arange(t_lo, t_hi, step)
    ok = np.array([ideal_peak_accel(r, v, T, g) <= a_max for T in Ts])
    if not ok.any():
        raise ValueError("no feasible time-to-go in range")
    i = np.argmax(ok)
    if i == 0:
        return Ts[0]
    lo, hi = Ts[i - 1], Ts[i]
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if ideal_peak_accel(r, v, mid, g) <= a_max:
            hi = mid
        else:
            lo = mid
    return hi
