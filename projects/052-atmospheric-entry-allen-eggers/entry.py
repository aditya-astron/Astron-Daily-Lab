# Astron Built this AP
# Project 052 | Thermal analysis / Atmospheric entry | Ballistic and lifting Earth entry: Allen-Eggers, equilibrium glide and Sutton-Graves heating
# Author: Aditya (aditya-astron)
# Date: 2026-11-21 IST
"""
Planar point-mass atmospheric entry over a spherical, non-rotating Earth with
an isothermal exponential atmosphere, plus the classical closed forms it is
checked against.

State y = [V, gamma, h, s, Q]: speed, flight-path angle (negative =
descending), altitude, downrange arc along the surface, and the stagnation
heat load integrated so far.

    dV/dt     = -D/m - g sin(gamma)
    dgamma/dt = L/(m V) - (g/V - V/r) cos(gamma)
    dh/dt     = V sin(gamma)
    ds/dt     = V cos(gamma) Re / r

    D/m = rho V^2 / (2 beta),   beta = m / (CD A)      ballistic coefficient
    L/m = (L/D) D/m
    rho = rho0 exp(-h/H),  g = mu / r^2,  r = Re + h

Stagnation-point convective heating is Sutton-Graves, q = k sqrt(rho/rn) V^3.

Closed forms (gravity and curvature dropped, gamma frozen at entry):

  Allen-Eggers ballistic, with s = sin|gamma_e|
    V(rho)      = Ve exp(-rho H / (2 beta s))
    peak decel  = Ve^2 s / (2 e H)  at rho = beta s / H,      V = Ve e^(-1/2)
    peak q      = k Ve^3 sqrt(beta s / (3 e H rn)) at rho = beta s / (3H), V = Ve e^(-1/6)
    heat load   = k Ve^2 sqrt(pi beta H / (rn s)) erf(sqrt(rho_g H / (beta s)))

  Equilibrium glide (gamma ~ 0, lift balances weight less centrifugal relief)
    L/m = g (1 - V^2/Vc^2),   Vc^2 = g r
    decel = g (1 - V^2/Vc^2) / (L/D)
    range = (r L/D / 2) ln[(1 - Vf^2/Vc^2) / (1 - V0^2/Vc^2)]
    sin(gamma_eq) = -2 H g / (V^2 L/D)
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
from scipy.special import erf

MU_E = 3.986004418e14   # m^3/s^2
R_E = 6.371e6           # m, mean radius
G0 = 9.80665            # m/s^2
RHO0 = 1.225            # kg/m^3
H_SCALE = 7200.0        # m, single-scale-height fit over 0-80 km
K_SG = 1.7415e-4        # Sutton-Graves constant for air, kg^0.5 / m


@dataclass(frozen=True)
class Vehicle:
    beta: float          # ballistic coefficient m/(CD A), kg/m^2
    lift_to_drag: float  # L/D, 0 for ballistic
    nose_radius: float   # m


@dataclass(frozen=True)
class Planet:
    mu: float = MU_E
    radius: float = R_E
    rho0: float = RHO0
    scale_height: float = H_SCALE
    gravity: bool = True     # switch off to recover the Allen-Eggers model
    curvature: bool = True   # switch off for a flat-Earth straight-line path

    def rho(self, h):
        return self.rho0 * np.exp(-np.asarray(h) / self.scale_height)

    def g(self, h):
        if not self.gravity:
            return 0.0
        return self.mu / (self.radius + h) ** 2


@dataclass
class Trajectory:
    veh: Vehicle
    pl: Planet
    sol: object          # solve_ivp result with dense output
    t: np.ndarray
    V: np.ndarray
    gamma: np.ndarray
    h: np.ndarray
    s: np.ndarray
    decel: np.ndarray    # aerodynamic load magnitude sqrt(D^2 + L^2)/m, m/s^2
    qdot: np.ndarray     # stagnation heat flux, W/m^2
    skipped: bool = False  # left the atmosphere again instead of being captured

    def loads(self, t):
        V, _, h, _, _ = self.sol.sol(t)
        rho = self.pl.rho(h)
        d = rho * V * V / (2 * self.veh.beta) * np.hypot(1.0, self.veh.lift_to_drag)
        q = K_SG * np.sqrt(rho / self.veh.nose_radius) * V ** 3
        return d, q

    def peak(self, name):
        """Peak value refined on the dense interpolant, with h and V there."""
        which = {"decel": 0, "qdot": 1}[name]
        i = int(np.argmax(getattr(self, name)))
        lo, hi = self.t[max(i - 1, 0)], self.t[min(i + 1, len(self.t) - 1)]
        res = minimize_scalar(lambda t: -self.loads(t)[which], bounds=(lo, hi),
                              method="bounded", options={"xatol": 1e-9})
        V, _, h, _, _ = self.sol.sol(res.x)
        return -res.fun, h, V

    def heat_load(self):
        """Integrated stagnation heat load, J/m^2 (carried as a fifth ODE state)."""
        return float(self.sol.y[4, -1])


def rhs(t, y, veh: Vehicle, pl: Planet):
    V, gam, h = y[0], y[1], y[2]
    r = pl.radius + h
    rho = pl.rho(h)
    d = rho * V * V / (2.0 * veh.beta)
    lift = veh.lift_to_drag * d
    g = pl.g(h)
    cg, sg = np.cos(gam), np.sin(gam)
    centrifugal = V / r if pl.curvature else 0.0
    dV = -d - g * sg
    dgam = lift / V - (g / V - centrifugal) * cg
    dh = V * sg
    ds = V * cg * (pl.radius / r if pl.curvature else 1.0)
    dQ = K_SG * np.sqrt(rho / veh.nose_radius) * V ** 3
    return [dV, dgam, dh, ds, dQ]


def fly(veh: Vehicle, pl: Planet, V0, gamma0, h0=120e3, h_end=0.0,
        V_end=150.0, rtol=1e-11, atol=1e-9, n_out=8001):
    """Integrate until the vehicle reaches h_end or slows to V_end."""
    def ground(t, y, *a):
        return y[2] - h_end
    ground.terminal, ground.direction = True, -1

    def slow(t, y, *a):
        return y[0] - V_end
    slow.terminal, slow.direction = True, -1

    def skip(t, y, *a):          # climbed back out through the interface
        return y[2] - (h0 + 1e3)
    skip.terminal, skip.direction = True, 1

    sol = solve_ivp(rhs, (0.0, 1e5), [V0, gamma0, h0, 0.0, 0.0], args=(veh, pl),
                    method="DOP853", rtol=rtol, atol=atol,
                    events=(ground, slow, skip), dense_output=True)
    if sol.status < 0:
        raise RuntimeError(sol.message)
    t = np.linspace(0.0, sol.t[-1], n_out)
    V, gam, h, s, _ = sol.sol(t)
    tr = Trajectory(veh, pl, sol, t, V, gam, h, s, None, None,
                    skipped=sol.t_events[2].size > 0)
    tr.decel, tr.qdot = tr.loads(t)
    return tr


# ---------------------------------------------------------------- closed forms

def ae_velocity(rho, Ve, beta, gamma_e, rho_e=0.0, H=H_SCALE):
    s = np.sin(abs(gamma_e))
    return Ve * np.exp(-(rho - rho_e) * H / (2.0 * beta * s))


def ae_peak_decel(Ve, beta, gamma_e, H=H_SCALE, rho0=RHO0):
    """Peak deceleration (m/s^2) and the altitude where it occurs."""
    s = np.sin(abs(gamma_e))
    a = Ve * Ve * s / (2.0 * np.e * H)
    h = -H * np.log(beta * s / (H * rho0))
    return a, h


def ae_peak_heating(Ve, beta, gamma_e, rn, H=H_SCALE, rho0=RHO0, k=K_SG):
    """Peak Sutton-Graves flux (W/m^2), its altitude and the speed there."""
    s = np.sin(abs(gamma_e))
    rho = beta * s / (3.0 * H)
    V = Ve * np.exp(-1.0 / 6.0)
    return k * np.sqrt(rho / rn) * V ** 3, -H * np.log(rho / rho0), V


def ae_heat_load(Ve, beta, gamma_e, rn, rho_g=RHO0, rho_e=0.0, H=H_SCALE, k=K_SG):
    """Integrated stagnation heat load (J/m^2) between densities rho_e and rho_g.

    Ve is the speed at rho = 0. The integrand goes as rho^(-1/2), so the thin
    upper layer is not negligible: starting at 120 km instead of rho = 0
    drops about 2 sqrt(x_e / pi) of the load, 0.5 % for a 5 deg entry.
    """
    s = np.sin(abs(gamma_e))
    x = lambda rho: np.sqrt(rho * H / (beta * s))
    return (k * Ve * Ve * np.sqrt(np.pi * beta * H / (rn * s))
            * (erf(x(rho_g)) - erf(x(rho_e))))


def ae_speed_at_zero_density(V_start, beta, gamma_e, rho_start, H=H_SCALE):
    """Back out the rho = 0 speed the closed forms use from a finite start altitude."""
    s = np.sin(abs(gamma_e))
    return V_start * np.exp(rho_start * H / (2.0 * beta * s))


def eq_glide_gamma(V, lift_to_drag, h, pl: Planet = Planet()):
    g = pl.g(h)
    return np.arcsin(-2.0 * pl.scale_height * g / (V * V * lift_to_drag))


def eq_glide_altitude(V, veh: Vehicle, pl: Planet = Planet(), h_guess=60e3):
    """Altitude where L/m = g - V^2/r; fixed-point on r because g and r move with h."""
    h = h_guess
    beta_L = veh.beta / veh.lift_to_drag      # m / (CL A)
    for _ in range(50):
        r = pl.radius + h
        rho = 2.0 * beta_L * (pl.g(h) - V * V / r) / (V * V)
        h_new = -pl.scale_height * np.log(rho / pl.rho0)
        if abs(h_new - h) < 1e-6:
            break
        h = h_new
    return h


def eq_glide_range(V0, Vf, lift_to_drag, r, g):
    vc2 = g * r
    return 0.5 * r * lift_to_drag * np.log((1 - Vf * Vf / vc2) / (1 - V0 * V0 / vc2))
