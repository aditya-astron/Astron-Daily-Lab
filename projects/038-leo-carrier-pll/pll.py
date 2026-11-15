# Astron Built this AP
# Project 038 | Signal processing / RF | Carrier-tracking PLL for an S-band LEO pass
# Author: Aditya (aditya-astron)
# Date: 2026-11-15 IST
"""
Ground-receiver carrier PLL for a 2.25 GHz downlink from a 500 km circular LEO.

The loop runs on 1 ms integrate-and-dump correlations, so it is modelled at
the phase level: each update the phase detector sees

    z_k = exp(j (phi_k - phihat_k)) + n_k,   E|n_k|^2 = 1 / ((C/N0) T)
    e_k = angle(z_k)                          (four-quadrant, residual carrier)

and at high SNR var(e_k) = 1 / (2 (C/N0) T). The filter is built from
rectangular (forward Euler) accumulators, the way it would sit in an FPGA:

    order 1:  NCO rate = K e
    order 2:  w += T wn^2 e;            NCO rate = w + 2 zeta wn e
    order 3:  a += T w0^3 e;  w += T (a + c2 w0^2 e);  NCO rate = w + c1 w0 e

Every accumulator chain ends in w0^n, so the steady-state error to a constant
n-th phase derivative is exactly phi^(n) / w0^n in the discrete loop too.

Loop noise bandwidth (one-sided) of the discrete loop is taken from its own
impulse response h (noise in, NCO phase out): sigma^2 = var(e) sum(h^2) =
B_L / (C/N0)  gives  B_L = sum(h^2) / (2 T).

Pass geometry: non-rotating Earth, ground station in the orbit plane (an
overhead pass, which is the worst case for Doppler rate), elevation >= 10 deg.
"""

import numpy as np

MU = 3.986004418e14
RE = 6378137.0
C_LIGHT = 299792458.0

# Third-order loop filter F(s) = c1 w0 + c2 w0^2 / s + w0^3 / s^2 with the
# Kaplan & Hegarty coefficients (their b3 = 2.4, a3 = 1.1): B_L = 0.7845 w0.
C1, C2 = 2.4, 1.1
ZETA = 1.0 / np.sqrt(2.0)


# ---------------------------------------------------------------- loop design

def analog_bl(order, w0):
    """Closed-form one-sided noise bandwidth [Hz] of the continuous loop."""
    if order == 1:
        return w0 / 4.0
    if order == 2:
        return 0.5 * w0 * (ZETA + 1.0 / (4.0 * ZETA))
    return w0 * (C1**2 * C2 + C2**2 - C1) / (4.0 * (C1 * C2 - 1.0))


def w0_for_bl(order, bl):
    return bl / analog_bl(order, 1.0)


class Loop:
    """Vectorised digital PLL: one state vector entry per Monte Carlo path."""

    def __init__(self, order, bl, T, npaths=1):
        self.order, self.T = order, T
        self.w0 = w0_for_bl(order, bl)
        self.ph = np.zeros(npaths)   # NCO phase [rad]
        self.w = np.zeros(npaths)    # frequency accumulator [rad/s]
        self.a = np.zeros(npaths)    # frequency-rate accumulator [rad/s^2]

    def step(self, e):
        T, w0 = self.T, self.w0
        if self.order == 1:
            rate = w0 * e
        elif self.order == 2:
            self.w += T * w0**2 * e
            rate = self.w + 2.0 * ZETA * w0 * e
        else:
            self.a += T * w0**3 * e
            self.w += T * (self.a + C2 * w0**2 * e)
            rate = self.w + C1 * w0 * e
        self.ph += T * rate

    def transition_matrix(self):
        """Linear map of [ph, w, a] for a zero input phase (e = -ph)."""
        T, w0 = self.T, self.w0
        if self.order == 1:
            return np.array([[1.0 - T * w0]])
        if self.order == 2:
            k1, k2 = 2.0 * ZETA * w0, w0**2
            # w' = w - T k2 ph ; ph' = ph + T (w' - k1 ph)
            return np.array([[1.0 - T * k1 - T * T * k2, T],
                             [-T * k2, 1.0]])
        k1, k2, k3 = C1 * w0, C2 * w0**2, w0**3
        return np.array([
            [1.0 - T * k1 - T * T * k2 - T**3 * k3, T, T * T],
            [-T * k2 - T * T * k3, 1.0, T],
            [-T * k3, 0.0, 1.0]])


def impulse_response(order, bl, T, n=None):
    """NCO phase response to a unit impulse on the phase-detector output."""
    if n is None:
        n = int(60.0 / (bl * T)) + 200
    lp = Loop(order, bl, T)
    h = np.empty(n)
    for k in range(n):
        h[k] = lp.ph[0]
        lp.step(np.array([1.0 if k == 0 else 0.0]) - lp.ph)
    return h


def digital_bl(order, bl, T):
    """Actual B_L of the discrete loop: sum(h^2) / 2T via a Lyapunov solve.

    With e = v - ph the update is x' = A x + G v, and since ph enters only
    through e, G = e1 - A e1. sum h^2 = e1' P e1 with P = A P A' + G G'.
    """
    from scipy.linalg import solve_discrete_lyapunov
    A = Loop(order, bl, T).transition_matrix()
    G = -A[:, :1].copy()
    G[0, 0] += 1.0
    P = solve_discrete_lyapunov(A, G @ G.T)
    return P[0, 0] / (2.0 * T)


def spectral_radius(order, bl, T):
    return np.max(np.abs(np.linalg.eigvals(Loop(order, bl, T).transition_matrix())))


# ------------------------------------------------------------- pass geometry

class OverheadPass:
    def __init__(self, alt=500e3, fc=2.25e9, el_min_deg=10.0):
        self.r = RE + alt
        self.n = np.sqrt(MU / self.r**3)
        self.fc = fc
        self.k = 2.0 * np.pi * fc / C_LIGHT      # rad of carrier per metre
        # Half-pass angle where elevation hits el_min: r cos(th) - R = rho sin(el)
        el = np.radians(el_min_deg)
        th_e = np.arccos(RE * np.cos(el) / self.r) - el
        self.t_half = th_e / self.n

    def range_derivs(self, t):
        """rho and its first three time derivatives (closed form)."""
        R, r, n = RE, self.r, self.n
        th = n * t
        rho = np.sqrt(R * R + r * r - 2.0 * R * r * np.cos(th))
        r1 = R * r * n * np.sin(th) / rho
        r2 = (R * r * n * n * np.cos(th) - r1**2) / rho
        r3 = (-R * r * n**3 * np.sin(th) - 3.0 * r1 * r2) / rho
        return rho, r1, r2, r3

    def phase(self, t):
        """Received carrier phase [rad] and its derivatives, one-way downlink."""
        rho, r1, r2, r3 = self.range_derivs(t)
        rho0 = self.r - RE
        return -self.k * (rho - rho0), -self.k * r1, -self.k * r2, -self.k * r3

    def elevation_deg(self, t):
        rho = self.range_derivs(t)[0]
        return np.degrees(np.arcsin((self.r * np.cos(self.n * t) - RE) / rho))

    def max_derivs(self, npts=200001):
        t = np.linspace(-self.t_half, self.t_half, npts)
        _, p1, p2, p3 = self.phase(t)
        return np.abs(p1).max(), np.abs(p2).max(), np.abs(p3).max()


# -------------------------------------------------------------- simulation

def track(order, bl, T, t, phi, cn0_dbhz, npaths=1, seed=0, init=None):
    """Run the loop on a phase history phi(t). Returns error histories and slips.

    The tracking error is phi - phihat, unwrapped; a cycle slip is a change in
    round(err / 2 pi). init = (w0, a0) seeds the accumulators, standing in for
    the FLL/Doppler-predict acquisition that precedes PLL lock.
    """
    rng = np.random.default_rng(seed)
    lp = Loop(order, bl, T, npaths)
    lp.ph[:] = phi[0]
    if init is not None:
        lp.w[:], lp.a[:] = init
    cn0 = 10.0 ** (cn0_dbhz / 10.0)
    sig = np.sqrt(1.0 / (2.0 * cn0 * T)) if np.isfinite(cn0) else 0.0
    err = np.empty((len(t), npaths))
    cycle = np.zeros(npaths)
    slips = np.zeros(npaths, dtype=int)
    for k in range(len(t)):
        d = phi[k] - lp.ph
        err[k] = d
        c = np.round(d / (2.0 * np.pi))
        slips += (c != cycle)
        cycle = c
        if sig > 0.0:
            z = np.exp(1j * d) + sig * (rng.standard_normal(npaths)
                                        + 1j * rng.standard_normal(npaths))
            e = np.angle(z)
        else:
            e = np.angle(np.exp(1j * d))
        lp.step(e)
    return err, slips


# ----------------------------------------------------- first-order slip theory

def mean_slip_time_closed(rho, bl):
    """Viterbi/Lindsey: first-order PLL mean time to a cycle slip [s]."""
    from scipy.special import i0
    return np.pi**2 * rho * i0(rho) ** 2 / (2.0 * bl)


def mean_slip_time_quad(rho, bl, n=20001):
    """Same quantity from the Fokker-Planck first-passage double integral.

    dphi = -K sin(phi) dt + sqrt(2D) dW with D = K / rho; mean exit time from
    (-2 pi, 2 pi) starting at 0 equals reflect-at-0 / absorb-at-2 pi:
        T = (1/D) int_0^2pi e^{-rho cos y} int_0^y e^{rho cos z} dz dy
    """
    from scipy.integrate import cumulative_trapezoid, simpson
    K = 4.0 * bl
    D = K / rho
    y = np.linspace(0.0, 2.0 * np.pi, n)
    inner = cumulative_trapezoid(np.exp(rho * np.cos(y)), y, initial=0.0)
    return simpson(np.exp(-rho * np.cos(y)) * inner, x=y) / D


def simulate_first_order_slips(rho, bl, npaths, dt_k=0.004, seed=1, tmax_mult=12.0):
    """Euler-Maruyama of the continuous first-order loop until |phi| >= 2 pi."""
    rng = np.random.default_rng(seed)
    K = 4.0 * bl
    dt = dt_k / K
    s = np.sqrt(2.0 * (K / rho) * dt)
    phi = np.zeros(npaths)
    t_exit = np.full(npaths, np.nan)
    alive = np.arange(npaths)
    nmax = int(tmax_mult * mean_slip_time_closed(rho, bl) / dt)
    for k in range(1, nmax + 1):
        p = phi[alive]
        p += -K * np.sin(p) * dt + s * rng.standard_normal(p.size)
        phi[alive] = p
        out = np.abs(p) >= 2.0 * np.pi
        if out.any():
            t_exit[alive[out]] = k * dt
            alive = alive[~out]
            if alive.size == 0:
                break
    return t_exit


def stationary_cos_mean(rho, bl, npaths=4000, tsim_bl=40.0, dt_k=0.004, seed=2):
    """E[cos phi] of the first-order loop after burn-in; Tikhonov gives I1/I0."""
    rng = np.random.default_rng(seed)
    K = 4.0 * bl
    dt = dt_k / K
    s = np.sqrt(2.0 * (K / rho) * dt)
    phi = np.zeros(npaths)
    nsteps = int(tsim_bl / bl / dt)
    acc, cnt = 0.0, 0
    for k in range(nsteps):
        phi += -K * np.sin(phi) * dt + s * rng.standard_normal(npaths)
        if k > nsteps // 4 and k % 50 == 0:
            acc += np.cos(phi).mean()
            cnt += 1
    return acc / cnt


# ----------------------------------------------------------------- design

def budget_deg(order, bl, cn0_dbhz, dyn, T=None):
    """Rule-of-thumb tracking budget 3 sigma_thermal + e_dyn [deg], limit 45.

    With T given the thermal term uses the discrete loop's actual B_L, and an
    unstable discrete loop scores infinity.
    """
    cn0 = 10.0 ** (cn0_dbhz / 10.0)
    b_eff = bl
    if T is not None:
        if spectral_radius(order, bl, T) >= 1.0:
            return np.inf
        b_eff = digital_bl(order, bl, T)
    e_dyn = dyn / w0_for_bl(order, bl) ** order
    return np.degrees(3.0 * np.sqrt(b_eff / cn0) + e_dyn)


def optimal_bl(order, cn0_dbhz, dyn, T=None):
    """B_L minimising the budget.

    Continuous loop: d/dB [3 sqrt(B/C) + g B^-n] = 0 gives
    B^(n+1/2) = (2 n / 3) g sqrt(C), g = dyn (B_L/w0)^n. Discrete loop: numeric,
    bracketed around the continuous optimum.
    """
    cn0 = 10.0 ** (cn0_dbhz / 10.0)
    g = dyn * analog_bl(order, 1.0) ** order
    b0 = (2.0 * order / 3.0 * g * np.sqrt(cn0)) ** (1.0 / (order + 0.5))
    if T is None:
        return b0
    from scipy.optimize import minimize_scalar
    f = lambda lb: budget_deg(order, np.exp(lb), cn0_dbhz, dyn, T)
    hi = np.log(min(4.0 * b0, 0.25 / T))
    lo = min(np.log(0.25 * b0), hi - np.log(16.0))
    r = minimize_scalar(f, bounds=(lo, hi), method="bounded",
                        options={"xatol": 1e-7})
    return np.exp(r.x)


def threshold_cn0(order, dyn, T=None, limit_deg=45.0):
    """Lowest C/N0 at which the optimally tuned loop meets the budget."""
    from scipy.optimize import brentq
    f = lambda c: budget_deg(order, optimal_bl(order, c, dyn, T), c, dyn, T) - limit_deg
    return brentq(f, 5.0, 80.0, xtol=1e-6)


def pass_monte_carlo(order, bl, T, cn0_dbhz, npaths, ps=None, seed=7):
    """Full pass from 10 deg to 10 deg. Returns (paths with a slip, noise sigma).

    The noise sigma is the rms of (err - noise-free err) over slip-free paths,
    to compare with the linear sqrt(B_L / (C/N0)).
    """
    ps = ps or OverheadPass()
    t = np.arange(-ps.t_half, ps.t_half, T)
    phi, ph1, ph2, _ = ps.phase(t)
    init = (ph1[0], ph2[0])
    e0, _ = track(order, bl, T, t, phi, np.inf, init=init)
    err, slips = track(order, bl, T, t, phi, cn0_dbhz, npaths, seed, init=init)
    ok = slips == 0
    sig = np.sqrt(np.mean((err[:, ok] - e0) ** 2)) if ok.any() else np.nan
    return int(np.count_nonzero(slips)), sig, np.abs(e0).max()


def main():
    ps = OverheadPass()
    p1, p2, p3 = ps.max_derivs()
    print("S-band 2.25 GHz, 500 km circular orbit, overhead pass above 10 deg")
    print(f"  pass length           {2 * ps.t_half:8.1f} s")
    print(f"  max Doppler           {p1 / 2 / np.pi / 1e3:8.2f} kHz")
    print(f"  max Doppler rate      {p2 / 2 / np.pi:8.1f} Hz/s   (at TCA)")
    print(f"  max Doppler jerk      {p3 / 2 / np.pi:8.2f} Hz/s^2")

    print("\nDigital loop noise bandwidth vs design (T = 1 ms)")
    print("  order  B_L design  B_L actual  B_L*T  |pole|max")
    for order in (2, 3):
        for bl in (5.0, 20.0, 50.0, 100.0):
            print(f"  {order:5d}  {bl:9.1f}  {digital_bl(order, bl, 1e-3):10.2f}"
                  f"  {bl * 1e-3:6.3f}  {spectral_radius(order, bl, 1e-3):8.4f}")

    print("\nDesign: 3 sigma + e_dyn <= 45 deg at the optimal B_L")
    print("  order  T [ms]  C/N0 thr  B_L [Hz]  B_L act  e_dyn [deg]  CN0*T [dB]")
    cases = []
    for order, dyn, T in ((2, p2, 1e-3), (3, p3, 1e-3), (3, p3, 1e-2)):
        c = threshold_cn0(order, dyn, T)
        b = optimal_bl(order, c, dyn, T)
        cases.append((order, T, c, b))
        print(f"  {order:5d}  {T * 1e3:6.0f}  {c:8.2f}  {b:8.2f}  {digital_bl(order, b, T):7.2f}"
              f"  {np.degrees(dyn / w0_for_bl(order, b) ** order):11.2f}"
              f"  {c + 10 * np.log10(T):10.2f}")

    n = 40
    print(f"\nFull-pass Monte Carlo at the design B_L ({n} passes per point)")
    print("  order  T [ms]  C/N0  passes w/ slip  sigma sim  sigma lin  max|e| no noise")
    for order, T, c, b in cases:
        lin = lambda cc: np.degrees(np.sqrt(digital_bl(order, b, T) / 10 ** (cc / 10)))
        for dc in (-2.0, 0.0, 2.0):
            k, sig, emax = pass_monte_carlo(order, b, T, c + dc, n)
            s = f"{np.degrees(sig):9.2f}" if np.isfinite(sig) else f"{'-':>9}"
            print(f"  {order:5d}  {T * 1e3:6.0f}  {c + dc:5.1f}  {k:9d}/{n}  "
                  f"  {s}  {lin(c + dc):9.2f}  {np.degrees(emax):9.2f}")


if __name__ == "__main__":
    main()
