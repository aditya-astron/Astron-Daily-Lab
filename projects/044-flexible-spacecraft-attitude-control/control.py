# Astron Built this AP
# Project 044 | PID and control systems | Attitude control of a spacecraft with flexible solar arrays
# Author: Aditya (aditya-astron)
# Date: 2026-11-17 IST
"""
Sampled-data PD attitude loop with an optional structural filter.

Timing: hub angle and gyro rate are sampled at t_k, the command is computed
and held (ZOH) from t_{k+1}: one full sample of transport delay, which is the
usual assumption for a flight computer that reads sensors at the start of a
minor frame and writes the wheel command at the start of the next.

Rigid-body gains from a continuous second-order target (wn, zeta):
    Kp = J_t wn^2,   Kd = 2 zeta wn J_t
Loop transfer, broken at the actuator (negative feedback, char. eq. 1 + L):
    L(z) = F(z) z^-1 [Kp P_theta(z) + Kd P_rate(z)]
"""

from dataclasses import dataclass

import numpy as np
from scipy.linalg import expm


def c2d_zoh(A, B, T):
    n, m = B.shape
    E = expm(np.block([[A, B], [np.zeros((m, n + m))]]) * T)
    return E[:n, :n], E[:n, n:]


@dataclass
class Filter:
    """Discrete SISO filter in state-space form (A, B, C, D)."""
    A: np.ndarray
    B: np.ndarray
    C: np.ndarray
    D: float
    name: str = "none"

    def freq(self, z):
        n = self.A.shape[0]
        if n == 0:
            return np.full_like(z, self.D, dtype=complex)
        out = np.empty(z.shape, dtype=complex)
        for i, zi in enumerate(z.ravel()):
            out.ravel()[i] = (self.C @ np.linalg.solve(zi * np.eye(n) - self.A, self.B)).item() + self.D
        return out


def unity():
    return Filter(np.zeros((0, 0)), np.zeros((0, 1)), np.zeros((1, 0)), 1.0, "none")


def tf2ss_biquad(b, a, name):
    """Direct-form biquad b(z)/a(z) (a[0] = 1) to controllable canonical SS."""
    b = np.asarray(b, float) / a[0]
    a = np.asarray(a, float) / a[0]
    A = np.array([[-a[1], -a[2]], [1.0, 0.0]])
    B = np.array([[1.0], [0.0]])
    C = np.array([[b[1] - b[0] * a[1], b[2] - b[0] * a[2]]])
    return Filter(A, B, C, b[0], name)


def bilinear_biquad(num, den, T, w_warp):
    """Bilinear transform of (n0 s^2 + n1 s + n2)/(d0 s^2 + d1 s + d2),
    prewarped so the response at w_warp is exact."""
    c = w_warp / np.tan(w_warp * T / 2.0)
    def poly(p):
        p0, p1, p2 = p
        # s = c (z - 1)/(z + 1); multiply through by (z + 1)^2, coefficients in z^-1
        return np.array([p0 * c**2 + p1 * c + p2,
                         -2 * p0 * c**2 + 2 * p2,
                         p0 * c**2 - p1 * c + p2])
    return poly(num), poly(den)


def notch(w0, T, zeta=0.5):
    """Full-depth notch (s^2 + w0^2)/(s^2 + 2 zeta w0 s + w0^2), width set by zeta.

    Below w0 it lags (up to -90 deg), above w0 it leads: that sign change is
    what decides which side of the notch a structural mode can sit on.
    """
    b, a = bilinear_biquad([1, 0, w0**2], [1, 2 * zeta * w0, w0**2], T, w0)
    return tf2ss_biquad(b, a, f"notch {w0:.2f} rad/s, zeta {zeta}")


def lowpass2(wf, T, zeta=0.7):
    b, a = bilinear_biquad([0, 0, wf**2], [1, 2 * zeta * wf, wf**2], T, wf)
    return tf2ss_biquad(b, a, f"low-pass {wf:.2f} rad/s, zeta {zeta}")


class Loop:
    def __init__(self, sat, T, wn, zeta, filt=None, n_modes=2, freq_scale=1.0, gain=1.0):
        self.T = T
        self.Kp = sat.J_t * wn**2 * gain
        self.Kd = 2 * zeta * wn * sat.J_t * gain
        self.filt = filt if filt is not None else unity()
        A, B, C = sat.state_space(n_modes, freq_scale)
        self.A, self.B, self.C = A, B, C
        self.Ad, self.Bd = c2d_zoh(A, B, T)
        self.K = np.array([[self.Kp, self.Kd]])

    def closed_loop(self):
        """Augmented state [x, x_f, u_held]; returns (Phi, G_ref) with ref = theta_cmd."""
        Ad, Bd, f = self.Ad, self.Bd, self.filt
        nx, nf = Ad.shape[0], f.A.shape[0]
        KC = self.K @ self.C                       # e_k = -KC x + Kp ref
        n = nx + nf + 1
        Phi = np.zeros((n, n))
        G = np.zeros(n)
        Phi[:nx, :nx] = Ad
        Phi[:nx, -1:] = Bd                         # plant driven by the held command
        Phi[nx:nx + nf, :nx] = -f.B @ KC
        Phi[nx:nx + nf, nx:nx + nf] = f.A
        G[nx:nx + nf] = f.B[:, 0] * self.Kp
        Phi[-1, :nx] = -f.D * KC[0]                # command computed now, applied next frame
        Phi[-1, nx:nx + nf] = f.C[0]
        G[-1] = f.D * self.Kp
        return Phi, G

    def spectral_radius(self):
        return np.max(np.abs(np.linalg.eigvals(self.closed_loop()[0])))

    def loop_tf(self, w):
        z = np.exp(1j * w * self.T)
        n = self.Ad.shape[0]
        P = np.array([self.C @ np.linalg.solve(zi * np.eye(n) - self.Ad, self.Bd)[:, 0] for zi in z])
        return self.filt.freq(z) / z * (P @ self.K[0])

    def margins(self, n=200000):
        """Nyquist crossings on a log grid up to the Nyquist frequency.

        Returns (phase_cross, gain_cross): phase_cross lists (w, 1/|L|) where the
        loop phase passes -180 deg (mod 360), i.e. the loop-gain multipliers at
        which a closed-loop pole reaches the unit circle; gain_cross lists
        (w, phase margin deg) where |L| = 1.
        """
        w = np.logspace(-3, np.log10(np.pi / self.T * 0.9999), n)
        L = self.loop_tf(w)
        ph = np.unwrap(np.angle(L))
        mag = np.abs(L)
        k = np.floor((ph + np.pi) / (2 * np.pi))
        phase_cross = []
        for i in np.nonzero(np.diff(k) != 0)[0]:
            target = 2 * np.pi * max(k[i], k[i + 1]) - np.pi
            a = (target - ph[i]) / (ph[i + 1] - ph[i])
            phase_cross.append((w[i] + a * (w[i + 1] - w[i]),
                                1.0 / (mag[i] + a * (mag[i + 1] - mag[i]))))
        gain_cross = []
        for i in np.nonzero(np.diff(np.sign(mag - 1.0)) != 0)[0]:
            a = (1.0 - mag[i]) / (mag[i + 1] - mag[i])
            p = ph[i] + a * (ph[i + 1] - ph[i])
            gain_cross.append((w[i] + a * (w[i + 1] - w[i]),
                               np.degrees((p + np.pi) % (2 * np.pi))))
        return phase_cross, gain_cross

    def mode_table(self):
        """Closed-loop poles mapped back to s = ln(z)/T: (|s| rad/s, damping ratio)."""
        ev = np.linalg.eigvals(self.closed_loop()[0])
        ev = ev[(np.abs(ev) > 1e-9)]
        s = np.log(ev.astype(complex)) / self.T
        s = s[s.imag >= -1e-12]
        return sorted((abs(x), -x.real / abs(x)) for x in s)

    def simulate(self, ref, N):
        Phi, G = self.closed_loop()
        X = np.zeros((N + 1, Phi.shape[0]))
        for k in range(N):
            X[k + 1] = Phi @ X[k] + G * ref
        return X
