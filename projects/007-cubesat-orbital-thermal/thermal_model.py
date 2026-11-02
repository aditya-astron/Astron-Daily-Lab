# Astron Built this AP
# Project 007 | Thermal analysis | Transient orbital thermal model of a nadir-pointing 3U CubeSat
# Author: Aditya (aditya-astron)
# Date: 2026-11-02 IST
"""Seven-node lumped thermal model of a 3U CubeSat in a circular LEO.

Nodes: the six external panels (+X ram, -X wake, +Y, -Y, +Z nadir, -Z zenith)
and one internal node lumping the electronics stack. Each panel sees

    Q_abs = alpha A (S cos(theta_s) [sunlit]  +  a S F cos(psi))  +  eps A q_IR F
    Q_out = eps sigma A T^4                    (deep space at ~0 K)

where F is the panel-to-Earth view factor, psi the angle between the sun and
the sub-satellite point, a the Bond albedo and q_IR Earth's outgoing longwave
flux. Panels exchange heat with each other through edge joints (linear
conductances) and with the internal node through standoffs (conductance) and
the cavity (radiation, sigma * R_ij * (T_j^4 - T_i^4)).

The node equations C dT/dt = Q(t) + G T + radiation are stiff-ish and
nonlinear, so they are marched with backward Euler and a Newton solve at each
step. The run repeats whole orbits until the temperature history is periodic.
"""

from dataclasses import dataclass, field
import math

import numpy as np

SIGMA = 5.670374419e-8          # W m^-2 K^-4
MU_EARTH = 3.986004418e14       # m^3 s^-2
R_EARTH = 6371.0e3              # m, mean radius
SOLAR_CONST = 1361.0            # W/m^2 at 1 AU
ALBEDO = 0.30
EARTH_IR = 237.0                # W/m^2, global mean OLR

NODE_NAMES = ["+X", "-X", "+Y", "-Y", "+Z", "-Z", "internal"]
# Outward normals in the body frame. Body frame is LVLH-aligned:
# x = along velocity, z = nadir, y completes the right-handed set (anti orbit normal).
FACE_NORMALS = np.array([
    [1, 0, 0], [-1, 0, 0],
    [0, 1, 0], [0, -1, 0],
    [0, 0, 1], [0, 0, -1],
], dtype=float)


def earth_view_factor(lam, h_ratio):
    """View factor from a flat plate to a sphere.

    lam      angle between the plate normal and the nadir direction [rad]
    h_ratio  H = r / R_earth (>1)

    Closed form for a differential plate tilted from the local vertical
    (Howell's configuration factor catalog, case B-22; also Gilmore,
    Spacecraft Thermal Control Handbook, ch. 2). With sin(phi) = 1/H:
      lam <= pi/2 - phi : F = cos(lam) / H^2               (whole disk visible)
      lam >= pi/2 + phi : F = 0                            (Earth behind plate)
      otherwise         : partial-visibility expression below
    """
    H = h_ratio
    phi = math.asin(1.0 / H)
    if lam <= math.pi / 2 - phi:
        return math.cos(lam) / H**2
    if lam >= math.pi / 2 + phi:
        return 0.0
    x = math.sqrt(H * H - 1.0)
    s, c = math.sin(lam), math.cos(lam)
    term1 = 0.5 - math.asin(min(1.0, x / (H * s))) / math.pi
    inner = max(0.0, 1.0 - H * H * c * c)
    term2 = (c * math.acos(max(-1.0, min(1.0, -x * c / s))) - x * math.sqrt(inner)) / (math.pi * H * H)
    return term1 + term2


def in_cylindrical_shadow(r_hat, sun_hat, r):
    """True when the spacecraft is inside Earth's cylindrical shadow."""
    d = float(np.dot(r_hat, sun_hat))
    if d >= 0.0:
        return False
    return r * math.sqrt(max(0.0, 1.0 - d * d)) < R_EARTH


def eclipse_fraction_analytic(alt, beta):
    """Fraction of a circular orbit spent in the cylindrical shadow.

    f = acos( sqrt(h^2 + 2 R h) / (r cos beta) ) / pi  for |beta| < beta*,
    beta* = asin(R / r); zero above that.
    """
    r = R_EARTH + alt
    if abs(beta) >= math.asin(R_EARTH / r):
        return 0.0
    return math.acos(math.sqrt(alt * alt + 2 * R_EARTH * alt) / (r * math.cos(beta))) / math.pi


@dataclass
class Spacecraft:
    """Default numbers describe a 10 x 10 x 34 cm 3U with body-mounted cells on
    the +/-Y and -Z faces, bare anodised aluminium elsewhere."""
    area: np.ndarray = field(default_factory=lambda: np.array(
        [0.010, 0.010, 0.034, 0.034, 0.034, 0.034]))                  # m^2
    alpha: np.ndarray = field(default_factory=lambda: np.array(
        [0.25, 0.25, 0.80, 0.80, 0.25, 0.80]))                         # solar absorptance
    eps: np.ndarray = field(default_factory=lambda: np.array(
        [0.80, 0.80, 0.85, 0.85, 0.80, 0.85]))                         # IR emittance
    # Heat capacities: 1.6 mm Al panels (+ cells), and the internal stack.
    cap: np.ndarray = field(default_factory=lambda: np.array(
        [45.0, 45.0, 160.0, 160.0, 150.0, 160.0, 1800.0]))             # J/K
    g_edge: float = 0.25       # W/K, joint conductance between adjacent panels
    g_standoff: float = 0.12   # W/K, each panel to the internal stack
    eps_internal: float = 0.6  # effective cavity emittance (panels <-> stack)
    q_internal: float = 4.0    # W, average electronics dissipation

    def conductance_matrix(self):
        n = 7
        G = np.zeros((n, n))
        opposite = {0: 1, 1: 0, 2: 3, 3: 2, 4: 5, 5: 4}
        for i in range(6):
            for j in range(6):
                if i != j and opposite[i] != j:
                    G[i, j] = self.g_edge
            G[i, 6] = G[6, i] = self.g_standoff
        np.fill_diagonal(G, -G.sum(axis=1))
        return G

    def radiation_coupling(self):
        """sigma * R_ij, panel interior faces to the stack. Radiative exchange
        factor approximated as eps_eff * A_panel (stack fully enclosed)."""
        R = np.zeros((7, 7))
        for i in range(6):
            R[i, 6] = R[6, i] = SIGMA * self.eps_internal * self.area[i]
        return R


@dataclass
class Orbit:
    alt: float = 500e3                 # m
    beta: float = math.radians(20.0)   # sun angle above the orbit plane

    @property
    def radius(self):
        return R_EARTH + self.alt

    @property
    def period(self):
        return 2 * math.pi * math.sqrt(self.radius**3 / MU_EARTH)


def environment_loads(sc: Spacecraft, orb: Orbit, theta):
    """Absorbed external heat per panel at orbit angle theta (measured from the
    sub-solar point projected into the orbit plane). Returns (Q[7], sunlit)."""
    r = orb.radius
    sun = np.array([math.cos(orb.beta), 0.0, math.sin(orb.beta)])
    r_hat = np.array([math.cos(theta), math.sin(theta), 0.0])
    v_hat = np.array([-math.sin(theta), math.cos(theta), 0.0])
    z_b = -r_hat
    y_b = np.cross(z_b, v_hat)
    dcm = np.vstack([v_hat, y_b, z_b])      # inertial -> body
    sun_b = dcm @ sun

    sunlit = not in_cylindrical_shadow(r_hat, sun, r)
    cos_psi = float(np.dot(r_hat, sun))      # sun elevation at sub-satellite point
    h_ratio = r / R_EARTH

    Q = np.zeros(7)
    for i in range(6):
        n = FACE_NORMALS[i]
        lam = math.acos(max(-1.0, min(1.0, n[2])))   # nadir is +z in body frame
        F = earth_view_factor(lam, h_ratio)
        direct = SOLAR_CONST * max(0.0, float(n @ sun_b)) if sunlit else 0.0
        albedo = ALBEDO * SOLAR_CONST * F * max(0.0, cos_psi)
        Q[i] = sc.area[i] * (sc.alpha[i] * (direct + albedo) + sc.eps[i] * EARTH_IR * F)
    Q[6] = sc.q_internal
    return Q, sunlit


def backward_euler_step(T, dt, Q, cap, G, Rint, emit, tol=1e-9, max_iter=30):
    """Solve cap*(T1 - T)/dt = Q + G T1 + Rint-exchange(T1) - emit*T1^4 with Newton.

    emit[i] = eps_i sigma A_i for panels (radiation to space), 0 for internal.
    """
    T1 = T.copy()
    for _ in range(max_iter):
        T4 = T1**4
        rad = Rint @ T4 - Rint.sum(axis=1) * T4
        f = cap * (T1 - T) / dt - Q - G @ T1 - rad + emit * T4
        J = np.diag(cap / dt) - G
        J -= 4.0 * Rint * T1**3                                   # d/dTj of R_ij Tj^4
        J += np.diag(4.0 * (Rint.sum(axis=1) + emit) * T1**3)
        dT = np.linalg.solve(J, -f)
        T1 += dT
        if np.max(np.abs(dT)) < tol:
            return T1
    raise RuntimeError("Newton did not converge in backward Euler step")


def simulate(sc=None, orb=None, steps_per_orbit=720, max_orbits=60,
             periodic_tol=1e-3, T0=290.0):
    """March whole orbits until the start-of-orbit state repeats within
    periodic_tol [K]. Returns a dict with the final orbit's history."""
    sc = sc or Spacecraft()
    orb = orb or Orbit()
    P = orb.period
    dt = P / steps_per_orbit
    G = sc.conductance_matrix()
    Rint = sc.radiation_coupling()
    emit = np.zeros(7)
    emit[:6] = sc.eps * SIGMA * sc.area

    # Environment depends only on orbit angle, so tabulate it once.
    thetas = 2 * math.pi * np.arange(1, steps_per_orbit + 1) / steps_per_orbit
    loads = [environment_loads(sc, orb, th) for th in thetas]
    Qtab = np.array([q for q, _ in loads])
    sunlit = np.array([s for _, s in loads])

    T = np.full(7, T0)
    for orbit in range(1, max_orbits + 1):
        T_start = T.copy()
        hist = np.empty((steps_per_orbit, 7))
        q_in = q_out = 0.0
        for k in range(steps_per_orbit):
            T = backward_euler_step(T, dt, Qtab[k], sc.cap, G, Rint, emit)
            hist[k] = T
            q_in += Qtab[k].sum() * dt
            q_out += float(emit @ T**4) * dt
        if np.max(np.abs(T - T_start)) < periodic_tol:
            break
    else:
        raise RuntimeError("no periodic solution within max_orbits")

    return {
        "time": thetas / (2 * math.pi) * P,
        "T": hist,
        "sunlit": sunlit,
        "orbits": orbit,
        "period": P,
        "energy_in": q_in,
        "energy_out": q_out,
        "dt": dt,
    }


def report(res):
    T = res["T"] - 273.15
    ecl = 1.0 - res["sunlit"].mean()
    lines = [
        f"orbit period {res['period']:.1f} s, eclipse fraction {ecl:.3f}, "
        f"periodic after {res['orbits']} orbits",
        f"energy balance over final orbit: in {res['energy_in']/1e3:.2f} kJ, "
        f"out {res['energy_out']/1e3:.2f} kJ",
        "",
        f"{'node':>9} {'min C':>8} {'mean C':>8} {'max C':>8} {'swing K':>8}",
    ]
    for i, name in enumerate(NODE_NAMES):
        col = T[:, i]
        lines.append(f"{name:>9} {col.min():8.1f} {col.mean():8.1f} {col.max():8.1f} "
                     f"{col.max() - col.min():8.1f}")
    return "\n".join(lines)


if __name__ == "__main__":
    for beta_deg in (0.0, 20.0, 60.0, 75.0):
        print(f"--- 500 km circular, beta = {beta_deg:.0f} deg ---")
        print(report(simulate(orb=Orbit(beta=math.radians(beta_deg)))))
        print()
