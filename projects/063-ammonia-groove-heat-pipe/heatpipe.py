# Astron Built this AP
# Project 063 | Thermal analysis | Axially grooved aluminium/ammonia heat pipe transport limits
# Author: Aditya (aditya-astron)
# Date: 2026-11-25 IST
"""
Transport limits of an axially grooved aluminium/ammonia heat pipe.

Geometry: N rectangular grooves of width w and depth d cut into the wall
around a vapour core of diameter Dv (the fin-tip circle). Fins between grooves
are t wide at the tip, so N = floor(pi Dv / (w + t)).

Capillary limit (steady, the one that sizes groove pipes):

    2 sigma / w  >=  dP_liq + dP_vap + rho_l g (Lt sin(phi) + Dv cos(phi))

phi is the adverse tilt (evaporator up). The Dv cos(phi) term is the head
needed to keep the top grooves wet in 1 g and is dropped in 0 g. The flat
meniscus sits at the condenser end cap, the tightest one (radius w/2) at the
evaporator end cap.

Liquid in a groove with a shear-free free surface is half of a symmetric
w x 2d duct, so the exact Poiseuille series for a rectangle applies.
"""
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
import CoolProp.CoolProp as CP

G0 = 9.80665


# --------------------------------------------------------------- properties
@dataclass(frozen=True)
class Sat:
    T: float        # K
    P: float        # Pa
    rho_l: float    # kg/m^3
    rho_v: float
    mu_l: float     # Pa s
    mu_v: float
    sigma: float    # N/m
    hfg: float      # J/kg

    @property
    def merit(self):
        """Liquid transport factor sigma rho_l hfg / mu_l, W/m^2."""
        return self.sigma * self.rho_l * self.hfg / self.mu_l


@lru_cache(maxsize=512)
def ammonia(T):
    f = "Ammonia"
    p = lambda k, q: CP.PropsSI(k, "T", T, "Q", q, f)
    return Sat(T=T, P=p("P", 0), rho_l=p("D", 0), rho_v=p("D", 1),
               mu_l=p("V", 0), mu_v=p("V", 1), sigma=p("I", 0),
               hfg=p("H", 1) - p("H", 0))


# --------------------------------------------------------- groove hydraulics
def rect_flow_coeff(W, H, nterms=200):
    """Q = C * G / mu for fully developed flow in a W x H rectangle
    (G = -dP/dx). Series from White, Viscous Fluid Flow, eq. 3-48, with
    half-sides a = W/2 <= b = H/2 swapped as needed for fast convergence."""
    a, b = sorted((W / 2.0, H / 2.0))
    n = np.arange(1, 2 * nterms, 2, dtype=float)
    s = np.sum(np.tanh(n * np.pi * b / (2 * a)) / n**5)
    return 4.0 * b * a**3 / 3.0 * (1.0 - 192.0 * a / (np.pi**5 * b) * s)


def groove_flow_coeff(w, d):
    """Open rectangular groove, shear-free top: half of a w x 2d duct."""
    return 0.5 * rect_flow_coeff(w, 2.0 * d)


def shah_london_fre(alpha):
    """Fanning f*Re for a rectangular duct, aspect ratio alpha <= 1
    (Shah & London 1978 fit; within 0.06 % of the exact series)."""
    c = [1.0, -1.3553, 1.9467, -1.7012, 0.9564, -0.2537]
    return 24.0 * sum(ci * alpha**i for i, ci in enumerate(c))


def groove_flow_coeff_fd(w, d, n=120):
    """Finite-difference check of groove_flow_coeff: solve lap(u) = -1 on the
    w x d groove with no-slip on the walls and du/dz = 0 on the free surface,
    integrate u. Second-order, so two grids are Richardson-extrapolated."""
    from scipy.sparse import diags, kron, identity
    from scipy.sparse.linalg import spsolve

    def solve(nx):
        nz = max(4, int(round(nx * d / w)))
        hx, hz = w / nx, d / nz
        # interior x nodes 1..nx-1 (u=0 at the side walls)
        mx = nx - 1
        Dx = diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(mx, mx)) / hx**2
        # z nodes 1..nz; u=0 at z=0, mirror ghost at z=nz+1 for zero shear
        Dz = diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(nz, nz)).tolil()
        Dz[nz - 1, nz - 2] = 2.0
        Dz = Dz.tocsr() / hz**2
        A = kron(identity(nz), Dx) + kron(Dz, identity(mx))
        u = spsolve(A.tocsc(), -np.ones(mx * nz)).reshape(nz, mx)
        # trapezoid in z (top row weight 1/2, bottom wall u=0), x walls are 0
        wz = np.ones(nz); wz[-1] = 0.5
        return hx * hz * np.sum(u * wz[:, None])

    q1, q2 = solve(n), solve(2 * n)
    return (4 * q2 - q1) / 3


# ------------------------------------------------------------------ the pipe
@dataclass(frozen=True)
class Pipe:
    Dv: float = 8.5e-3     # vapour core / fin-tip diameter
    w: float = 0.5e-3      # groove width
    d: float = 1.2e-3      # groove depth
    t: float = 0.35e-3     # fin tip thickness (extrusion minimum)
    Le: float = 0.15       # evaporator
    La: float = 0.60       # adiabatic
    Lc: float = 0.30       # condenser

    @property
    def N(self):
        return int(np.floor(np.pi * self.Dv / (self.w + self.t)))

    @property
    def Lt(self):
        return self.Le + self.La + self.Lc

    @property
    def Leff(self):
        """Uniform evaporator/condenser: La + (Le + Lc)/2."""
        return self.La + 0.5 * (self.Le + self.Lc)

    @property
    def Av(self):
        return np.pi * self.Dv**2 / 4

    def liquid_resistance(self, s):
        """dP_liq per unit (Q * length), Pa/(W m)."""
        return s.mu_l / (s.rho_l * s.hfg * self.N * groove_flow_coeff(self.w, self.d))

    def vapour_resistance(self, s, Q=0.0):
        """Core pressure gradient per unit (Q * length), Pa/(W m). Laminar
        Hagen-Poiseuille (f Re = 16) below Re 2300, Blasius above; the peak
        Reynolds number at heat load Q decides. Near -40 C a 12.7 mm pipe
        sits right at the transition."""
        Re = 4.0 * Q / (s.hfg * np.pi * self.Dv * s.mu_v)
        fRe = 16.0 if Re < 2300.0 else 0.079 * Re**0.75
        return fRe / 16.0 * 128.0 * s.mu_v / (np.pi * s.rho_v * s.hfg * self.Dv**4)

    def cap_pressure(self, s):
        return 2.0 * s.sigma / self.w

    def gravity_head(self, s, tilt_mm=0.0, one_g=False):
        """Adverse elevation tilt_mm of the evaporator end over the condenser
        end; one_g adds the circumferential Dv cos(phi) ~ Dv term."""
        g = G0 if one_g else 0.0
        return s.rho_l * g * (tilt_mm * 1e-3 + self.Dv)

    def q_capillary(self, s, tilt_mm=0.0, one_g=False, Leff=None):
        Leff = self.Leff if Leff is None else Leff
        drive = self.cap_pressure(s) - self.gravity_head(s, tilt_mm, one_g)
        if drive <= 0.0:
            return 0.0
        rl = self.liquid_resistance(s)
        # Largest Q with Q (rl + rv(Q)) Leff <= drive. The left side rises
        # monotonically but jumps at the regime change, so: laminar answer if
        # it stays laminar, else the turbulent fixed point if it stays
        # turbulent, else the load that pins the core at Re = 2300.
        q_re = 2300.0 * s.hfg * np.pi * self.Dv * s.mu_v / 4.0
        q_lam = drive / ((rl + self.vapour_resistance(s, 0.0)) * Leff)
        if q_lam < q_re:
            return q_lam
        Q = q_lam
        for _ in range(100):
            Qn = drive / ((rl + self.vapour_resistance(s, max(Q, q_re))) * Leff)
            if abs(Qn - Q) < 1e-13 * Q:
                break
            Q = Qn
        return Qn if Qn >= q_re else q_re

    # Other classical limits (Faghri, Heat Pipe Science and Technology, ch. 4)
    def q_sonic(self, s):
        """Busse choked-vapour limit."""
        return 0.474 * self.Av * s.hfg * np.sqrt(s.rho_v * s.P)

    def q_viscous(self, s):
        """Busse vapour-pressure limit, relevant only near the triple point."""
        return self.Av * self.Dv**2 * s.hfg * s.rho_v * s.P / (64.0 * s.mu_v * self.Leff)

    def q_entrainment(self, s):
        """Shear-entrainment (Weber) limit, wick surface hydraulic radius w/2."""
        return self.Av * s.hfg * np.sqrt(s.sigma * s.rho_v / (2.0 * (self.w / 2)))

    def limits(self, T):
        s = ammonia(T)
        return {"capillary_0g": self.q_capillary(s),
                "sonic": self.q_sonic(s),
                "entrainment": self.q_entrainment(s),
                "viscous": self.q_viscous(s)}

    def vapour_state(self, s, Q):
        """Peak core Reynolds and Mach number at heat load Q."""
        m = Q / s.hfg
        Re = 4 * m / (np.pi * self.Dv * s.mu_v)
        c = CP.PropsSI("A", "T", s.T, "Q", 1, "Ammonia")
        Ma = m / (s.rho_v * self.Av * c)
        return Re, Ma


# ------------------------------------------------------- axial pressure march
def march(pipe, s, q_lin, Q, n=4001):
    """Integrate the liquid and vapour pressure gradients along the pipe for
    an arbitrary heat input/output distribution q_lin(x) (W/m, positive in,
    must integrate to zero). Returns x and the local capillary pressure
    P_v - P_l, referenced to a flat meniscus at the condenser end cap.

    Local axial heat flow Qx(x) = int_0^x q_lin; mass flow Qx/hfg runs as
    vapour towards the condenser and as liquid back, so both gradients add
    to the meniscus pressure difference."""
    from scipy.integrate import solve_ivp

    R = pipe.liquid_resistance(s) + pipe.vapour_resistance(s, Q)

    def rhs(x, y):
        Qx, dPc = y
        return [q_lin(x), -R * Qx]

    knots = [0, pipe.Le, pipe.Le + pipe.La, pipe.Lt]
    xs, ys = [], []
    y0 = [0.0, 0.0]
    for a, b in zip(knots[:-1], knots[1:]):
        sol = solve_ivp(rhs, (a, b), y0, rtol=1e-11, atol=1e-12,
                        dense_output=True, method="DOP853")
        xx = np.linspace(a, b, n // 3)
        xs.append(xx); ys.append(sol.sol(xx))
        y0 = sol.y[:, -1]
    x = np.concatenate(xs)
    y = np.concatenate(ys, axis=1)
    # shift so the condenser end cap has a flat meniscus
    dPc = y[1] - y[1][-1]
    return x, y[0], dPc


def uniform_profile(pipe, Q):
    qe, qc = Q / pipe.Le, Q / pipe.Lc

    def q(x):
        if x < pipe.Le:
            return qe
        if x < pipe.Le + pipe.La:
            return 0.0
        return -qc
    return q
