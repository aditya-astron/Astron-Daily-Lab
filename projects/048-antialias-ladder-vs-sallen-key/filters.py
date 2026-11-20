# Astron Built this AP
# Project 048 | Electronics and circuit analysis | ADC anti-alias filter: Chebyshev LC ladder vs. Sallen-Key cascade
# Author: Aditya (aditya-astron)
# Date: 2026-11-20 IST
"""
Building blocks for a 50 ohm anti-alias low-pass in front of a telemetry ADC.

Two realisations of the same Chebyshev type I response are compared:

* a doubly terminated LC ladder (shunt-C first, minimum inductor count),
  synthesised from the Matthaei closed-form prototype g-values;
* a cascade of unity-gain Sallen-Key biquads plus one passive RC pole,
  placed directly on the prototype poles.

Everything is evaluated on the j-omega axis from component values, so the
same code computes nominal responses, dissipation and tolerance Monte Carlo.
Arrays broadcast: component arrays may carry a leading trial axis.
"""

import math

import numpy as np

TWO_PI = 2.0 * math.pi


# ---------------------------------------------------------------------------
# Order selection and prototype
# ---------------------------------------------------------------------------

def _ripple_eps2(ap_db):
    return 10.0 ** (ap_db / 10.0) - 1.0


def cheby_order(ap_db, as_db, ws_over_wp):
    """Minimum Chebyshev I order for ripple ap_db and attenuation as_db at ws."""
    ratio = math.sqrt((10.0 ** (as_db / 10.0) - 1.0) / _ripple_eps2(ap_db))
    return math.ceil(math.acosh(ratio) / math.acosh(ws_over_wp))


def butter_order(ap_db, as_db, ws_over_wp):
    """Minimum Butterworth order with the same ap_db held at the passband edge."""
    ratio = math.sqrt((10.0 ** (as_db / 10.0) - 1.0) / _ripple_eps2(ap_db))
    return math.ceil(math.log(ratio) / math.log(ws_over_wp))


def cheby_g(n, ap_db):
    """
    Low-pass prototype element values g0..g_{n+1} (Matthaei, Young & Jones 4.05).

    beta  = ln coth(Ap / 17.37),  gamma = sinh(beta / 2n)
    a_k   = sin((2k-1) pi / 2n),  b_k   = gamma^2 + sin^2(k pi / n)
    g1    = 2 a1 / gamma,         g_k   = 4 a_{k-1} a_k / (b_{k-1} g_{k-1})
    g_{n+1} = 1 (n odd) or coth^2(beta/4) (n even).
    Normalised to a 1 ohm source and a ripple edge at 1 rad/s.
    """
    beta = math.log(1.0 / math.tanh(ap_db / (40.0 / math.log(10.0))))
    gamma = math.sinh(beta / (2 * n))
    a = [math.sin((2 * k - 1) * math.pi / (2 * n)) for k in range(1, n + 1)]
    b = [gamma ** 2 + math.sin(k * math.pi / n) ** 2 for k in range(1, n + 1)]
    g = [1.0, 2.0 * a[0] / gamma]
    for k in range(2, n + 1):
        g.append(4.0 * a[k - 2] * a[k - 1] / (b[k - 2] * g[k - 1]))
    g.append(1.0 if n % 2 else 1.0 / math.tanh(beta / 4.0) ** 2)
    return g


def butter_g(n):
    """Butterworth prototype g-values, 3 dB edge at 1 rad/s, equal terminations."""
    return [1.0] + [2.0 * math.sin((2 * k - 1) * math.pi / (2 * n))
                    for k in range(1, n + 1)] + [1.0]


def cheby_poles(n, ap_db):
    """Left-half-plane poles of the Chebyshev I prototype (ripple edge 1 rad/s)."""
    eps = math.sqrt(_ripple_eps2(ap_db))
    mu = math.asinh(1.0 / eps) / n
    k = np.arange(1, n + 1)
    theta = (2 * k - 1) * math.pi / (2 * n)
    return -math.sinh(mu) * np.sin(theta) + 1j * math.cosh(mu) * np.cos(theta)


def cheby_mag2(w_over_wp, n, ap_db):
    """|H|^2 = 1 / (1 + eps^2 T_n^2(w)), Chebyshev polynomial by cosh/cos form."""
    x = np.abs(np.asarray(w_over_wp, dtype=float))
    tn = np.where(x <= 1.0,
                  np.cos(n * np.arccos(np.minimum(x, 1.0))),
                  np.cosh(n * np.arccosh(np.maximum(x, 1.0))))
    return 1.0 / (1.0 + _ripple_eps2(ap_db) * tn ** 2)


def pole_group_delay(w, poles):
    """Group delay of an all-pole H(s) = prod(-p)/(s-p): sum of -Re p / |jw - p|^2."""
    w = np.asarray(w, dtype=float)[..., None]
    p = np.asarray(poles)
    return np.sum(-p.real / ((w - p.imag) ** 2 + p.real ** 2), axis=-1)


# ---------------------------------------------------------------------------
# Doubly terminated LC ladder
# ---------------------------------------------------------------------------

def ladder_values(g, r0, f_edge):
    """
    Denormalise a shunt-C-first prototype to r0 and a band edge f_edge (Hz).
    Odd positions (g1, g3, ...) are shunt capacitors, even ones series inductors.
    g_{n+1} is a load resistance after a shunt C but a load conductance after
    a series L (Matthaei's convention), which matters for even n.
    """
    w = TWO_PI * f_edge
    elems = []
    for k, gk in enumerate(g[1:-1], start=1):
        if k % 2:
            elems.append(("C", gk / (r0 * w)))
        else:
            elems.append(("L", gk * r0 / w))
    rl = r0 * g[-1] if elems[-1][0] == "C" else r0 / g[-1]
    return elems, rl


def ladder_sparams(f, kinds, values, rs, rl, q_l=math.inf, q_c=math.inf):
    """
    S21 and S11 of the ladder between rs and rl by ABCD cascade.

    Constant-Q loss model: an inductor is jwL + wL/Q_L, a capacitor
    jwC + wC/Q_C, the usual datasheet description near the operating band.
    values has shape (..., n_elements); f has shape (F,). Returns (..., F).
    S21 = 2 sqrt(rs/rl) V_L/V_s, so a lossless matched ladder has |S21| = 1.
    """
    w = TWO_PI * np.asarray(f, dtype=float)
    vals = np.asarray(values, dtype=float)[..., None, :]
    A = np.ones(vals.shape[:-1], dtype=complex)
    B = np.zeros_like(A)
    C = np.zeros_like(A)
    D = np.ones_like(A)
    for i, kind in enumerate(kinds):
        x = vals[..., i]
        if kind == "C":
            y = (1j * w + w / q_c) * x
            A, B, C, D = A + B * y, B, C + D * y, D
        else:
            z = (1j * w + w / q_l) * x
            A, B, C, D = A, A * z + B, C, C * z + D
    vl_over_vs = rl / (A * rl + B + rs * (C * rl + D))
    s21 = 2.0 * math.sqrt(rs / rl) * vl_over_vs
    zin = (A * rl + B) / (C * rl + D)
    s11 = (zin - rs) / (zin + rs)
    return s21, s11


# ---------------------------------------------------------------------------
# Sallen-Key cascade
# ---------------------------------------------------------------------------

def sallen_key_design(poles_rad_s, r):
    """
    Unity-gain Sallen-Key sections, equal resistors r, on given LHP poles.

    H = 1 / (s^2 R^2 C1 C2 + 2 s R C2 + 1)  =>  w0 = 1/(R sqrt(C1 C2)),
    Q = sqrt(C1/C2)/2, so C2 = 1/(2 Q w0 R) and C1 = 2 Q / (w0 R).
    C1 is the feedback capacitor, C2 the grounded one. A real pole becomes a
    passive R-C on the buffered output. Sections are ordered by rising Q,
    which keeps internal peaking (and op-amp swing) lowest.
    Returns a list of ("SK", R1, R2, C1, C2) and ("RC", R, C) tuples.
    """
    secs = []
    upper = [p for p in poles_rad_s if p.imag > 1e-9 * abs(p)]
    real = [p for p in poles_rad_s if abs(p.imag) <= 1e-9 * abs(p)]
    for p in sorted(upper, key=lambda p: abs(p) / (-2 * p.real)):
        w0 = abs(p)
        q = w0 / (-2.0 * p.real)
        secs.append(("SK", r, r, 2.0 * q / (w0 * r), 1.0 / (2.0 * q * w0 * r)))
    for p in real:
        secs.append(("RC", r, 1.0 / (-p.real * r)))
    return secs


def sallen_key_response(f, secs, gbw_hz=math.inf):
    """
    Cascade response. Each Sallen-Key stage is solved from its node equations
    with a single-pole op-amp A(s) = wt/s as follower, so finite gain-bandwidth
    is included exactly:
      node x: (Vi - Vx)/R1 + (Vp - Vx)/R2 + s C1 (Vo - Vx) = 0
      node p: (Vx - Vp)/R2 = s C2 Vp,  Vo = Vp A/(1 + A)
    Section parameters may be arrays with a leading trial axis.
    """
    s = 1j * TWO_PI * np.asarray(f, dtype=float)
    h = 1.0 + 0j
    for sec in secs:
        if sec[0] == "SK":
            r1, r2, c1, c2 = (np.asarray(v, dtype=float)[..., None] for v in sec[1:])
            if math.isinf(gbw_hz):
                beta = 1.0
            else:
                wt = TWO_PI * gbw_hz
                beta = wt / (s + wt)      # A/(1+A) with A = wt/s, safe at DC
            # Vp = Vx / (1 + s R2 C2); substitute into node x and solve Vx/Vi.
            kp = 1.0 / (1.0 + s * r2 * c2)
            g1, g2 = 1.0 / r1, 1.0 / r2
            vx = g1 / (g1 + g2 + s * c1 - g2 * kp - s * c1 * beta * kp)
            h = h * vx * kp * beta
        else:
            r, c = (np.asarray(v, dtype=float)[..., None] for v in sec[1:])
            h = h / (1.0 + s * r * c)
    return h
