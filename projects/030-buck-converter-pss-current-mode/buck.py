# Astron Built this AP
# Project 030 | Electronics and circuit analysis | Buck converter: exact periodic steady state, DCM boundary and peak-current-mode subharmonic stability
# Author: Aditya (aditya-astron)
# Date: 2026-11-11 IST
"""
Point-of-load buck converter off a 28 V spacecraft bus.

The switched circuit is linear between switching instants, so each interval
has an exact discrete solution x(t0+tau) = Phi(tau) x(t0) + Gamma(tau) Vin
from one augmented matrix exponential. Chaining the ON and OFF intervals gives
the one-period map, and its fixed point is the periodic steady state (PSS)
solved directly, with no start-up transient to wait out.

State x = [iL, vC]. Output node: vo = R (vC + rc iL) / (R + rc).

Three studies:
  1. Synchronous CCM buck: PSS vs. the textbook averaged results
     (Vo, inductor ripple, output ripple, efficiency).
  2. Non-synchronous buck at light load: DCM conversion ratio
     M = 2 / (1 + sqrt(1 + 4K/D^2)), K = 2L/(R T), and the CCM/DCM boundary
     K_crit = 1 - D.
  3. Peak-current-mode inner loop: cycle-to-cycle map of the inductor current,
     perturbation ratio alpha = -(m2 - ma)/(m1 + ma), subharmonic oscillation
     for D > 0.5 without slope compensation, and the boundary against
     Ridley's sampling-gain Q = 1 / (pi (mc D' - 0.5)).
"""
from dataclasses import dataclass, replace

import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq, fsolve


@dataclass(frozen=True)
class Buck:
    vin: float = 28.0      # bus voltage [V]
    L: float = 10e-6       # inductance [H]
    C: float = 47e-6       # output capacitance [F]
    R: float = 1.0         # load [ohm]
    rL: float = 15e-3      # inductor DCR [ohm]
    rc: float = 10e-3      # capacitor ESR [ohm]
    ron: float = 20e-3     # high- and low-side switch on-resistance [ohm]
    fs: float = 500e3      # switching frequency [Hz]

    @property
    def T(self):
        return 1.0 / self.fs

    def vo(self, x):
        x = np.asarray(x)
        return self.R * (x[..., 1] + self.rc * x[..., 0]) / (self.R + self.rc)

    def matrices(self, conducting=True):
        """A, B for an interval. conducting=False is the DCM idle interval
        (both switches off, iL held at zero, C discharges into R + rc)."""
        R, rc, L, C = self.R, self.rc, self.L, self.C
        p = R / (R + rc)
        if not conducting:
            A = np.array([[0.0, 0.0], [0.0, -1.0 / (C * (R + rc))]])
            return A, np.zeros(2)
        # L diL/dt = u Vin - (ron + rL) iL - vo,  C dvC/dt = (R iL - vC)/(R + rc)
        A = np.array([
            [-(self.ron + self.rL + rc * p) / L, -p / L],
            [R / (C * (R + rc)), -1.0 / (C * (R + rc))],
        ])
        B = np.array([1.0 / L, 0.0])
        return A, B


def step(A, B, u, tau):
    """Exact zero-order-hold solution over tau: returns (Phi, Gamma*u)."""
    n = A.shape[0]
    M = np.zeros((n + 1, n + 1))
    M[:n, :n] = A
    M[:n, n] = B * u
    E = expm(M * tau)
    return E[:n, :n], E[:n, n]


def propagate(A, B, u, x0, tau):
    Phi, g = step(A, B, u, tau)
    return Phi @ x0 + g


# ---------------------------------------------------------------- CCM, sync

def pss_sync(cv: Buck, D):
    """Fixed point of the synchronous-buck period map: x = Phi x + Gamma."""
    A, B = cv.matrices()
    P1, g1 = step(A, B, cv.vin, D * cv.T)
    P2, g2 = step(A, B, 0.0, (1 - D) * cv.T)
    Phi = P2 @ P1
    Gam = P2 @ g1 + g2
    return np.linalg.solve(np.eye(2) - Phi, Gam)


def waveform_sync(cv: Buck, D, n=4000):
    """Sample one PSS period on a uniform grid (exact at every sample)."""
    A, B = cv.matrices()
    x0 = pss_sync(cv, D)
    t = np.arange(n) * cv.T / n
    ton = D * cv.T
    xon_end = propagate(A, B, cv.vin, x0, ton)
    X = np.empty((n, 2))
    for k, tk in enumerate(t):
        if tk < ton:
            X[k] = propagate(A, B, cv.vin, x0, tk)
        else:
            X[k] = propagate(A, B, 0.0, xon_end, tk - ton)
    return t, X, ton


def interval_mean(f, A, B, u, x0, tau, nodes=24):
    """Exact-to-quadrature mean of f(x(t)) over an interval. The trajectory
    is a sum of decaying exponentials, so Gauss-Legendre converges fast."""
    g, w = np.polynomial.legendre.leggauss(nodes)
    s = 0.5 * tau * (g + 1.0)
    vals = [f(propagate(A, B, u, x0, sk)) for sk in s]
    return float(np.dot(w, vals) / 2.0)


def ccm_report(cv: Buck, D):
    A, B = cv.matrices()
    x0 = pss_sync(cv, D)
    x1 = propagate(A, B, cv.vin, x0, D * cv.T)

    def mean(f):
        return (D * interval_mean(f, A, B, cv.vin, x0, D * cv.T)
                + (1 - D) * interval_mean(f, A, B, 0.0, x1, (1 - D) * cv.T))

    _, X, _ = waveform_sync(cv, D)
    iL, vo = X[:, 0], cv.vo(X)
    pin = cv.vin * D * interval_mean(lambda x: x[0], A, B, cv.vin, x0, D * cv.T)
    pout = mean(lambda x: cv.vo(x) ** 2) / cv.R
    il_rms2 = mean(lambda x: x[0] ** 2)
    ic_rms2 = mean(lambda x: (x[0] - cv.vo(x) / cv.R) ** 2)
    return dict(
        vo_avg=mean(cv.vo), il_avg=mean(lambda x: x[0]),
        di_pp=iL.max() - iL.min(), dv_pp=vo.max() - vo.min(),
        il_rms2=il_rms2, ic_rms2=ic_rms2, pin=pin, pout=pout, eff=pout / pin,
        loss_ohmic=(cv.ron + cv.rL) * il_rms2 + cv.rc * ic_rms2,
    )


def ccm_analytic(cv: Buck, D):
    """Averaged-model predictions (ripple neglected in the DC terms)."""
    rs = cv.ron + cv.rL               # both switches have the same ron
    vo = D * cv.vin * cv.R / (cv.R + rs)
    io = vo / cv.R
    di = (cv.vin - vo - rs * io) * D * cv.T / cv.L
    dv_c = di / (8 * cv.fs * cv.C)    # capacitor charge ripple
    dv_esr = di * cv.rc
    eff = cv.R / (cv.R + rs)
    return dict(vo=vo, io=io, di=di, dv_c=dv_c, dv_esr=dv_esr, eff=eff)


# ------------------------------------------------------- non-sync, DCM/CCM

def period_map_async(cv: Buck, D, x0):
    """One period of a non-synchronous buck with an ideal diode.
    Returns (x_end, t_idle): t_idle > 0 means the inductor ran dry (DCM)."""
    A, B = cv.matrices(True)
    Az, Bz = cv.matrices(False)
    x1 = propagate(A, B, cv.vin, x0, D * cv.T)
    toff = (1 - D) * cv.T
    xe = propagate(A, B, 0.0, x1, toff)
    if xe[0] >= 0.0:
        return xe, 0.0
    # iL crosses zero inside the off interval: find the instant exactly.
    t0 = brentq(lambda s: propagate(A, B, 0.0, x1, s)[0], 0.0, toff, xtol=1e-15)
    xz = propagate(A, B, 0.0, x1, t0)
    xz[0] = 0.0
    return propagate(Az, Bz, 0.0, xz, toff - t0), toff - t0


def pss_async(cv: Buck, D):
    guess = np.array([0.0, dcm_ratio(cv, D) * cv.vin])
    x = fsolve(lambda x: period_map_async(cv, D, x)[0] - x, guess, xtol=1e-13)
    return x, period_map_async(cv, D, x)[1]


def dcm_ratio(cv: Buck, D):
    K = 2 * cv.L / (cv.R * cv.T)
    if K >= 1 - D:
        return D
    return 2.0 / (1.0 + np.sqrt(1.0 + 4.0 * K / D ** 2))


def vo_avg_async(cv: Buck, D, x0, n=4000):
    """Average output over one period by dense exact sampling."""
    A, B = cv.matrices(True)
    Az, _ = cv.matrices(False)
    x1 = propagate(A, B, cv.vin, x0, D * cv.T)
    toff = (1 - D) * cv.T
    _, tidle = period_map_async(cv, D, x0)
    tz = toff - tidle
    xz = propagate(A, B, 0.0, x1, tz)
    xz[0] = 0.0
    vs = []
    for tk in (np.arange(n) + 0.5) * cv.T / n:
        if tk < D * cv.T:
            x = propagate(A, B, cv.vin, x0, tk)
        elif tk < D * cv.T + tz:
            x = propagate(A, B, 0.0, x1, tk - D * cv.T)
        else:
            x = propagate(Az, B, 0.0, xz, tk - D * cv.T - tz)
        vs.append(cv.vo(x))
    return float(np.mean(vs))


# --------------------------------------------------- peak current mode

def pcm_map(i0, ic, m1, m2, ma, T):
    """Inductor current at the next clock edge for a peak-current modulator
    with a stiff output. The switch turns on at the clock and off when
    iL = ic - ma t. On-time is clamped to [0, T] (max duty)."""
    ton = np.clip((ic - i0) / (m1 + ma), 0.0, T)
    return i0 + m1 * ton - m2 * (T - ton)


def pcm_alpha(m1, m2, ma):
    return -(m2 - ma) / (m1 + ma)


def pcm_run(i0, ic, m1, m2, ma, T, n):
    seq = [i0]
    for _ in range(n):
        seq.append(pcm_map(seq[-1], ic, m1, m2, ma, T))
    return np.array(seq)


def pcm_fixed_point(ic, m1, m2, ma, T):
    """Valley current of the period-1 orbit: ton = D T with D = m2/(m1+m2)."""
    D = m2 / (m1 + m2)
    return ic - (m1 + ma) * D * T


def pcm_lyapunov(i0, ic, m1, m2, ma, T, n=20000, burn=500):
    """Largest Lyapunov exponent of the 1-D map. The map is piecewise linear:
    slope alpha where the comparator trips, slope 1 where the duty clamps at
    0 or 1, so lambda = (fraction of unclamped cycles) * ln|alpha|."""
    i, acc, used = i0, 0.0, 0
    a = abs(pcm_alpha(m1, m2, ma))
    for k in range(n + burn):
        ton = (ic - i) / (m1 + ma)
        if k >= burn:
            acc += np.log(a) if 0.0 < ton < T else 0.0
            used += 1
        i = pcm_map(i, ic, m1, m2, ma, T)
    return acc / used


def ridley_q(m1, m2, ma):
    D = m2 / (m1 + m2)
    mc = 1.0 + ma / m1
    den = np.pi * (mc * (1 - D) - 0.5)
    return np.inf if abs(den) < 1e-12 else 1.0 / den


def main():
    np.set_printoptions(precision=5)
    cv = Buck()
    D = 5.0 / 28.0 * (cv.R + cv.ron + cv.rL) / cv.R   # trim for 5.000 V
    sim, ana = ccm_report(cv, D), ccm_analytic(cv, D)
    print("1) Synchronous buck 28 V -> 5 V, 5 A, 500 kHz, exact PSS")
    print(f"   D = {D:.5f}")
    print(f"   Vo avg     PSS {sim['vo_avg']:.5f} V   averaged {ana['vo']:.5f} V")
    print(f"   dIL p-p    PSS {sim['di_pp']:.5f} A   analytic {ana['di']:.5f} A")
    print(f"   dVo p-p    PSS {sim['dv_pp']*1e3:.3f} mV  "
          f"cap {ana['dv_c']*1e3:.3f} mV + ESR {ana['dv_esr']*1e3:.3f} mV")
    print(f"   efficiency PSS {sim['eff']*100:.4f} %   averaged {ana['eff']*100:.4f} %")
    print(f"   energy balance  Pin - Pout = {(sim['pin']-sim['pout'])*1e3:.3f} mW,"
          f"  sum I^2 R = {sim['loss_ohmic']*1e3:.3f} mW")
    noesr = replace(cv, rc=0.0)
    s0, a0 = ccm_report(noesr, D), ccm_analytic(noesr, D)
    print(f"   rc = 0:  dVo p-p PSS {s0['dv_pp']*1e3:.4f} mV   dI/(8 fs C) {a0['dv_c']*1e3:.4f} mV")

    print("\n2) Non-synchronous lossless buck, light load, D = 0.25")
    lossless = replace(cv, rL=0.0, rc=0.0, ron=0.0)
    Dd = 0.25
    Rcrit = 2 * cv.L * cv.fs / (1 - Dd)
    print(f"   K_crit = 1 - D -> R_crit = {Rcrit:.3f} ohm")
    for R in (5.0, 20.0, 50.0, 200.0):
        c = replace(lossless, R=R)
        x0, tidle = pss_async(c, Dd)
        v = vo_avg_async(c, Dd, x0)
        print(f"   R = {R:6.1f}  mode {'DCM' if tidle > 0 else 'CCM'}  "
              f"M sim {v/c.vin:.5f}  analytic {dcm_ratio(c, Dd):.5f}")

    print("\n3) Peak current mode, 5 V -> 3.3 V POL (D = 0.66), 1 MHz")
    vin, vo, L, T = 5.0, 3.3, 2.2e-6, 1e-6
    m1, m2 = (vin - vo) / L, vo / L
    ic = 3.0
    for label, ma in (("no slope comp", 0.0),
                      ("ma = (m2-m1)/2 boundary", 0.5 * (m2 - m1)),
                      ("ma = m2/2", 0.5 * m2),
                      ("ma = m2 deadbeat", m2)):
        a = pcm_alpha(m1, m2, ma)
        istar = pcm_fixed_point(ic, m1, m2, ma, T)
        seq = pcm_run(istar + 1e-6, ic, m1, m2, ma, T, 4)
        meas = (seq[1] - istar) / (seq[0] - istar)
        q = ridley_q(m1, m2, ma)
        print(f"   {label:26s} alpha {a:+.4f} (map {meas:+.4f})  Ridley Q {q:+.3g}")
    i0 = pcm_fixed_point(ic, m1, m2, 0.0, T) + 1e-3
    seq = pcm_run(i0, ic, m1, m2, 0.0, T, 200)
    lam = pcm_lyapunov(i0, ic, m1, m2, 0.0, T)
    print(f"   no slope comp, cycles 195-200 valley currents: "
          + " ".join(f"{v:.3f}" for v in seq[-6:]) + " A")
    print(f"   Lyapunov exponent {lam:+.3f} per cycle (> 0: bounded, aperiodic, chaotic)")
    print(f"   with ma = m2/2: {pcm_lyapunov(i0, ic, m1, m2, 0.5*m2, T):+.3f} per cycle")


if __name__ == "__main__":
    main()
