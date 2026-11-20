# Astron Built this AP
# Project 050 | Rocket propulsion | Regenerative cooling of a LOX/methane thrust chamber
# Author: Aditya (aditya-astron)
# Date: 2026-11-20 IST
"""
Steady regenerative cooling of a milled-channel copper thrust chamber.

Hot gas side: Bartz correlation with the property-variation factor sigma,
local Mach from the isentropic area ratio and a turbulent recovery
temperature T_aw = T0 (1 + r (g-1)/2 M^2) / (1 + (g-1)/2 M^2), r = Pr^(1/3).

Wall: 1D conduction through the hot-wall ligament tw. The lands between
channels are straight fins of thickness b and height h, adiabatic at the
closeout, so each channel rejects

    q_c' = h_c (T_wc - T_b) (w + 2 eta h),   eta = tanh(m h)/(m h),
    m = sqrt(2 h_c / (k b))

per unit channel length while the gas deposits q_g * 2 pi r / N on it.

Coolant side: Gnielinski Nusselt number with a Sieder-Tate (mu_b/mu_w)^0.14
wall correction, which raises h_c while the bulk is liquid-like and lowers it
once the bulk is gas-like, Haaland friction with wall roughness, and
real-fluid methane properties from CoolProp (HEOS) on (h, p). Heat-transfer
deterioration near the pseudo-critical point and channel curvature effects
are not modelled.

The coolant enters at the nozzle exit and flows toward the injector
(counterflow). The march is Heun predictor-corrector in axial position with
the arc-length factor ds/dx = sqrt(1 + r'^2) applied to both heat input and
friction length.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq

R_UNIV = 8314.462618  # J/(kmol K)


# ---------------------------------------------------------------------------
# Isentropic flow
# ---------------------------------------------------------------------------
def area_ratio(M: float, g: float) -> float:
    """A/A* for Mach M."""
    t = (2.0 / (g + 1.0)) * (1.0 + 0.5 * (g - 1.0) * M * M)
    return t ** ((g + 1.0) / (2.0 * (g - 1.0))) / M


def mach_from_area(eps: float, g: float, supersonic: bool) -> float:
    if eps < 1.0 - 1e-12:
        raise ValueError("area ratio below 1")
    if abs(eps - 1.0) < 1e-12:
        return 1.0
    f = lambda M: area_ratio(M, g) - eps
    if supersonic:
        return brentq(f, 1.0 + 1e-12, 50.0, xtol=1e-14, rtol=1e-14)
    return brentq(f, 1e-9, 1.0 - 1e-12, xtol=1e-14, rtol=1e-14)


# ---------------------------------------------------------------------------
# Combustion gas (frozen, calorically perfect at chamber conditions)
# ---------------------------------------------------------------------------
@dataclass
class CombustionGas:
    T0: float          # chamber stagnation temperature, K
    gamma: float
    molar_mass: float  # kg/kmol
    pc: float          # chamber pressure, Pa

    @property
    def R(self) -> float:
        return R_UNIV / self.molar_mass

    @property
    def cp(self) -> float:
        return self.gamma * self.R / (self.gamma - 1.0)

    @property
    def cstar(self) -> float:
        g = self.gamma
        gam = (2.0 / (g + 1.0)) ** ((g + 1.0) / (2.0 * (g - 1.0)))
        return math.sqrt(self.R * self.T0 / g) / gam

    @property
    def prandtl(self) -> float:
        # Eucken-type estimate used by Bartz for combustion products
        return 4.0 * self.gamma / (9.0 * self.gamma - 5.0)

    @property
    def mu0(self) -> float:
        # Bartz: mu [lb/(in s)] = 46.6e-10 M^0.5 T[R]^0.6, converted to Pa s
        lb_in_s = 46.6e-10 * self.molar_mass ** 0.5 * (self.T0 * 1.8) ** 0.6
        return lb_in_s * 0.45359237 / 0.0254

    def recovery_temperature(self, M: float) -> float:
        g = self.gamma
        r = self.prandtl ** (1.0 / 3.0)
        s = 0.5 * (g - 1.0) * M * M
        return self.T0 * (1.0 + r * s) / (1.0 + s)


def bartz_sigma(Twg: float, T0: float, M: float, g: float) -> float:
    s = 1.0 + 0.5 * (g - 1.0) * M * M
    return 1.0 / ((0.5 * Twg / T0 * s + 0.5) ** 0.68 * s ** 0.12)


def bartz_htc(gas: CombustionGas, Dt: float, r_curv: float, area_ratio_local: float,
              M: float, Twg: float) -> float:
    """Bartz gas-side film coefficient, SI, W/(m^2 K)."""
    lead = (0.026 / Dt ** 0.2) * (gas.mu0 ** 0.2 * gas.cp / gas.prandtl ** 0.6) \
        * (gas.pc / gas.cstar) ** 0.8 * (Dt / r_curv) ** 0.1
    return lead * area_ratio_local ** -0.9 * bartz_sigma(Twg, gas.T0, M, gas.gamma)


# ---------------------------------------------------------------------------
# Chamber contour, throat at x = 0
# ---------------------------------------------------------------------------
@dataclass
class ChamberContour:
    """Cylinder, fillet, 30 deg cone, throat arcs 1.5 Rt / 0.382 Rt, 15 deg cone.

    The cylinder length is set by L* = V_chamber / A_t, where V_chamber is
    everything upstream of the throat.
    """
    Rt: float
    contraction: float = 4.0
    lstar: float = 1.0
    eps_exit: float = 8.0
    theta_conv: float = math.radians(30.0)
    theta_div: float = math.radians(15.0)
    ru_factor: float = 1.5
    rd_factor: float = 0.382
    fillet_factor: float = 1.0   # fillet radius / chamber radius

    def __post_init__(self):
        Rt, tc = self.Rt, self.theta_conv
        self.Rc = Rt * math.sqrt(self.contraction)
        self.Re = Rt * math.sqrt(self.eps_exit)
        self.Ru = self.ru_factor * Rt
        self.Rd = self.rd_factor * Rt
        self.R2 = self.fillet_factor * self.Rc * 0.5
        self.x1 = -self.Ru * math.sin(tc)
        self.r1 = Rt + self.Ru * (1.0 - math.cos(tc))
        r2 = self.Rc - self.R2 * (1.0 - math.cos(tc))
        if r2 <= self.r1:
            raise ValueError("contraction too small for the chosen arcs")
        self.x2 = self.x1 - (r2 - self.r1) / math.tan(tc)
        self.xf = self.x2 - self.R2 * math.sin(tc)
        # divergent arc end and cone
        td = self.theta_div
        self.x3 = self.Rd * math.sin(td)
        self.r3 = Rt + self.Rd * (1.0 - math.cos(td))
        self.x_exit = self.x3 + (self.Re - self.r3) / math.tan(td)
        # cylinder length from L*
        v_conv = self._volume(self.xf, 0.0)
        v_need = self.lstar * math.pi * Rt * Rt
        if v_need <= v_conv:
            raise ValueError("L* too small: converging section alone exceeds it")
        self.Lcyl = (v_need - v_conv) / (math.pi * self.Rc ** 2)
        self.x_inj = self.xf - self.Lcyl

    def _volume(self, a: float, b: float, n: int = 4000) -> float:
        x = np.linspace(a, b, n + 1)
        r = np.array([self.r(xi) for xi in x])
        return float(np.trapezoid(math.pi * r * r, x))

    def r(self, x: float) -> float:
        Rt = self.Rt
        if x <= self.xf:
            return self.Rc
        if x <= self.x2:
            return self.Rc - self.R2 + math.sqrt(max(self.R2 ** 2 - (x - self.xf) ** 2, 0.0))
        if x <= self.x1:
            return self.r1 + (self.x1 - x) * math.tan(self.theta_conv)
        if x <= 0.0:
            return Rt + self.Ru - math.sqrt(max(self.Ru ** 2 - x * x, 0.0))
        if x <= self.x3:
            return Rt + self.Rd - math.sqrt(max(self.Rd ** 2 - x * x, 0.0))
        return self.r3 + (x - self.x3) * math.tan(self.theta_div)

    def drdx(self, x: float) -> float:
        if x <= self.xf:
            return 0.0
        if x <= self.x2:
            u = x - self.xf
            return -u / math.sqrt(max(self.R2 ** 2 - u * u, 1e-30))
        if x <= self.x1:
            return -math.tan(self.theta_conv)
        if x <= 0.0:
            return x / math.sqrt(max(self.Ru ** 2 - x * x, 1e-30))
        if x <= self.x3:
            return x / math.sqrt(max(self.Rd ** 2 - x * x, 1e-30))
        return math.tan(self.theta_div)

    @property
    def x_start(self) -> float:
        return self.x_inj

    @property
    def x_end(self) -> float:
        return self.x_exit

    @property
    def throat_curvature(self) -> float:
        return 0.5 * (self.Ru + self.Rd)


# ---------------------------------------------------------------------------
# Coolant property models
# ---------------------------------------------------------------------------
@dataclass
class FluidState:
    T: float
    rho: float
    mu: float
    k: float
    cp: float


class CoolPropFluid:
    """Real-fluid properties on (h, p) via CoolProp HEOS."""

    def __init__(self, name: str = "Methane"):
        import CoolProp.CoolProp as CP
        self._CP = CP
        self.name = name
        self.AS = CP.AbstractState("HEOS", name)

    def _update(self, pair, a, b):
        # A failed flash can leave stale guesses in the state object, so a
        # fresh one is built and the update retried once before giving up.
        try:
            self.AS.update(pair, a, b)
        except ValueError:
            self.AS = self._CP.AbstractState("HEOS", self.name)
            self.AS.update(pair, a, b)

    def h_from_Tp(self, T: float, p: float) -> float:
        self._update(self._CP.PT_INPUTS, p, T)
        return self.AS.hmass()

    def state(self, h: float, p: float) -> FluidState:
        self._update(self._CP.HmassP_INPUTS, h, p)
        AS = self.AS
        return FluidState(AS.T(), AS.rhomass(), AS.viscosity(), AS.conductivity(), AS.cpmass())

    def viscosity(self, T: float, p: float) -> float:
        """Viscosity at the wall. Above the EOS validity limit a dilute-gas
        T^0.7 law is continued from the limit instead of extrapolating."""
        Tmax = self.AS.Tmax()
        self._update(self._CP.PT_INPUTS, p, min(T, Tmax))
        mu = self.AS.viscosity()
        return mu * (T / Tmax) ** 0.7 if T > Tmax else mu


class ConstantFluid:
    """Constant-property fluid, h = cp T. Used for the analytic check."""

    def __init__(self, cp: float, rho: float, mu: float, k: float):
        self.cp, self.rho, self.mu, self.k = cp, rho, mu, k

    def h_from_Tp(self, T: float, p: float) -> float:
        return self.cp * T

    def state(self, h: float, p: float) -> FluidState:
        return FluidState(h / self.cp, self.rho, self.mu, self.k, self.cp)

    def viscosity(self, T: float, p: float) -> float:
        return self.mu


# ---------------------------------------------------------------------------
# Coolant-side correlations
# ---------------------------------------------------------------------------
def haaland_f(Re: float, rel_rough: float) -> float:
    """Darcy friction factor, Haaland explicit form of Colebrook."""
    x = -1.8 * math.log10((rel_rough / 3.7) ** 1.11 + 6.9 / Re)
    return 1.0 / (x * x)


def gnielinski_nu(Re: float, Pr: float, f: float) -> float:
    return (f / 8.0) * (Re - 1000.0) * Pr / (1.0 + 12.7 * math.sqrt(f / 8.0) * (Pr ** (2.0 / 3.0) - 1.0))


# ---------------------------------------------------------------------------
# Jacket design and solver
# ---------------------------------------------------------------------------
@dataclass
class Jacket:
    """Milled channel jacket.

    Either the land width is fixed and the channel width takes up the rest of
    the pitch (`width=None`), or the channel width is fixed and the land grows
    with radius. `height` is a constant or a callable of the gas-side radius.
    """
    n_channels: int
    height: object         # channel height h, m, float or h(r)
    land: float            # land width b, m (minimum land if width is fixed)
    t_wall: float          # hot-wall ligament tw, m
    k_wall: float          # wall conductivity, W/(m K)
    width: float | None = None
    roughness: float = 1.6e-6
    wall_exp: float = 0.14  # Sieder-Tate exponent on mu_b / mu_w

    def channel(self, r: float):
        """Return (width w, land b, height h, flow area, hydraulic diameter) at radius r."""
        pitch = 2.0 * math.pi * (r + self.t_wall) / self.n_channels
        if self.width is None:
            w, b = pitch - self.land, self.land
        else:
            w, b = self.width, pitch - self.width
        if w <= 0.0 or b < self.land - 1e-12:
            raise ValueError("channel does not fit the pitch: reduce N, land or width")
        h = self.height(r) if callable(self.height) else self.height
        return w, b, h, w * h, 2.0 * w * h / (w + h)


@dataclass
class StationResult:
    Twg: float
    Twc: float
    q: float        # gas-side heat flux, W/m^2
    hg: float
    hc: float
    eta: float
    Taw: float
    M: float


def coolant_htc(jk: Jacket, st: FluidState, mdot_ch: float, Dh: float, area: float,
                mu_wall: float):
    G = mdot_ch / area
    Re = G * Dh / st.mu
    Pr = st.mu * st.cp / st.k
    f = haaland_f(Re, jk.roughness / Dh)
    Nu = gnielinski_nu(Re, Pr, f) * (st.mu / mu_wall) ** jk.wall_exp
    return Nu * st.k / Dh, f, Re


def solve_station(jk: Jacket, r: float, hg_fun, Taw: float, st: FluidState, mdot_ch: float,
                  mu_wall=lambda T: None):
    """mu_wall(T) gives coolant viscosity at the wall; None means bulk value."""
    # Balance gas convection = wall conduction = channel convection (per channel).
    muw = lambda T: mu_wall(T) or st.mu
    w, b, h, area, Dh = jk.channel(r)
    tw, kw = jk.t_wall, jk.k_wall
    perim_hot = 2.0 * math.pi * r / jk.n_channels

    def parts(Twg):
        hg = hg_fun(Twg)
        q = hg * (Taw - Twg)
        Twc = Twg - q * tw / kw
        return hg, q, Twc

    def resid(Twg):
        hg, q, Twc = parts(Twg)
        if Twc <= st.T:
            return -q * perim_hot  # wall colder than coolant: Twg too low
        hc, _, _ = coolant_htc(jk, st, mdot_ch, Dh, area, muw(Twc))
        m = math.sqrt(2.0 * hc / (kw * b))
        eta = math.tanh(m * h) / (m * h)
        return hc * (Twc - st.T) * (w + 2.0 * eta * h) - q * perim_hot

    Twg = brentq(resid, st.T + 1e-6, Taw - 1e-6, xtol=1e-9, rtol=1e-13, maxiter=200)
    hg, q, Twc = parts(Twg)
    hc, _, _ = coolant_htc(jk, st, mdot_ch, Dh, area, muw(Twc))
    m = math.sqrt(2.0 * hc / (kw * b))
    eta = math.tanh(m * h) / (m * h)
    return StationResult(Twg, Twc, q, hg, hc, eta, Taw, float("nan"))


@dataclass
class CoolingResult:
    x: np.ndarray
    r: np.ndarray
    M: np.ndarray
    hg: np.ndarray
    hc: np.ndarray
    q: np.ndarray
    Twg: np.ndarray
    Twc: np.ndarray
    Tb: np.ndarray
    p: np.ndarray
    h: np.ndarray
    Q_total: float = 0.0
    extras: dict = field(default_factory=dict)

    @property
    def dp(self) -> float:
        return float(self.p[0] - self.p[-1])


class RegenSolver:
    """Counterflow march from nozzle exit (coolant inlet) to injector.

    `gas_side(x, Twg)` returns (hg, Taw, M) and `geometry` exposes r, drdx,
    x_start (injector) and x_end (exit). Swapping either lets the same march
    be checked against closed-form cases.
    """

    def __init__(self, geometry, jacket: Jacket, fluid, mdot: float, gas_side):
        self.g = geometry
        self.jk = jacket
        self.fluid = fluid
        self.mdot = mdot
        self.gas_side = gas_side

    def _rates(self, x: float, h: float, p: float):
        """d(h)/d(s_flow) and d(p)/d(s_flow) per unit axial length along the flow."""
        g, jk = self.g, self.jk
        r = g.r(x)
        ds = math.sqrt(1.0 + g.drdx(x) ** 2)
        st = self.fluid.state(h, p)
        mdot_ch = self.mdot / jk.n_channels
        hg_fun = lambda Twg: self.gas_side(x, Twg)[0]
        _, Taw, M = self.gas_side(x, 1000.0)
        muw = lambda T: self.fluid.viscosity(T, p)
        res = solve_station(jk, r, hg_fun, Taw, st, mdot_ch, muw)
        res.M = M
        _, _, _, area, Dh = jk.channel(r)
        _, f, _ = coolant_htc(jk, st, mdot_ch, Dh, area, muw(res.Twc))
        G = mdot_ch / area
        dh = res.q * 2.0 * math.pi * r * ds / self.mdot
        dp = -f / Dh * G * G / (2.0 * st.rho) * ds
        return dh, dp, res, st

    def run(self, T_in: float, p_in: float, n: int = 400) -> CoolingResult:
        g = self.g
        xs = np.linspace(g.x_end, g.x_start, n + 1)   # flow direction
        h = self.fluid.h_from_Tp(T_in, p_in)
        p = p_in
        out = {k: [] for k in ("x", "r", "M", "hg", "hc", "q", "Twg", "Twc", "Tb", "p", "h")}

        def record(x, res, st, h, p):
            for k, v in (("x", x), ("r", g.r(x)), ("M", res.M), ("hg", res.hg), ("hc", res.hc),
                         ("q", res.q), ("Twg", res.Twg), ("Twc", res.Twc), ("Tb", st.T),
                         ("p", p), ("h", h)):
                out[k].append(v)

        dh0, dp0, res0, st0 = self._rates(xs[0], h, p)
        record(xs[0], res0, st0, h, p)
        Q = 0.0
        for i in range(n):
            dx = abs(xs[i + 1] - xs[i])
            hp, pp = h + dh0 * dx, p + dp0 * dx
            if pp < 0.5 * p_in:
                raise ValueError("jacket pressure drop exceeds half the inlet pressure")
            dh1, dp1, res1, st1 = self._rates(xs[i + 1], hp, pp)
            h += 0.5 * (dh0 + dh1) * dx
            p += 0.5 * (dp0 + dp1) * dx
            Q += 0.5 * (dh0 + dh1) * dx * self.mdot
            dh0, dp0, res0, st0 = self._rates(xs[i + 1], h, p)
            record(xs[i + 1], res0, st0, h, p)
        arr = {k: np.asarray(v) for k, v in out.items()}
        return CoolingResult(Q_total=Q, **arr)


def bartz_gas_side(gas: CombustionGas, contour: ChamberContour):
    """Gas-side closure for a chamber contour: x -> (hg(Twg), Taw, M)."""
    Dt = 2.0 * contour.Rt
    rc = contour.throat_curvature
    cache = {}

    def fn(x: float, Twg: float):
        if x not in cache:
            eps = (contour.r(x) / contour.Rt) ** 2
            M = mach_from_area(max(eps, 1.0), gas.gamma, supersonic=x > 0.0)
            cache[x] = (eps, M, gas.recovery_temperature(M))
        eps, M, Taw = cache[x]
        return bartz_htc(gas, Dt, rc, eps, M, Twg), Taw, M

    return fn
