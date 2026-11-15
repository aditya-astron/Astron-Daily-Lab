# Astron Built this AP
# Project 038 | Signal processing / RF | Carrier-tracking PLL for an S-band LEO pass
# Author: Aditya (aditya-astron)
# Date: 2026-11-15 IST
"""Verification for pll.py. Run: python3 test_pll.py"""

import numpy as np
from scipy.integrate import quad
from scipy.special import i0, i1

from pll import (RE, C_LIGHT, C1, C2, ZETA, Loop, OverheadPass, analog_bl,
                 budget_deg, digital_bl, impulse_response, mean_slip_time_closed,
                 mean_slip_time_quad, optimal_bl, simulate_first_order_slips,
                 spectral_radius, stationary_cos_mean, track, w0_for_bl)


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name} {detail}")
    assert ok, name


def closed_loop_h(order, s):
    if order == 1:
        num, den = 1.0, s + 1.0
    elif order == 2:
        num = 2 * ZETA * s + 1.0
        den = s * s + num
    else:
        num = C1 * s * s + C2 * s + 1.0
        den = s**3 + num
    return num / den


def test_analog_bandwidth():
    # B_L = (1/2 pi) int_0^inf |H(jw)|^2 dw, w0 = 1
    for order, ref in ((1, 0.25), (2, 0.53033), (3, 0.7845)):
        bl = quad(lambda w: abs(closed_loop_h(order, 1j * w)) ** 2, 0, np.inf,
                  limit=400)[0] / (2 * np.pi)
        check(f"analog B_L order {order}", abs(bl - analog_bl(order, 1.0)) < 1e-8
              and abs(bl - ref) < 5e-5, f"{bl:.6f} w0")


def test_digital_bandwidth():
    for order in (1, 2, 3):
        h = impulse_response(order, 20.0, 1e-3)
        b_imp = np.sum(h**2) / 2e-3
        b_lyap = digital_bl(order, 20.0, 1e-3)
        check(f"Lyapunov B_L = sum(h^2)/2T, order {order}",
              abs(b_imp / b_lyap - 1) < 1e-9, f"{b_lyap:.4f} Hz")
        b_small = digital_bl(order, 1.0, 1e-4)      # B_L T = 1e-4
        check(f"digital -> analog as B_L T -> 0, order {order}",
              abs(b_small - 1.0) < 2e-3, f"{b_small:.5f} Hz")
    check("design points stable", spectral_radius(3, 6.5, 1e-2) < 1
          and spectral_radius(2, 90.0, 1e-3) < 1)
    # Stability edges. Order 1: pole 1 - K T, edge at K T = 2 (B_L T = 0.5).
    # Order 2, Jury: 1 + tr + det = 4 - 2 T k1 - T^2 k2 > 0 -> w0 T = sqrt6 - sqrt2.
    from scipy.optimize import brentq
    edge = {o: brentq(lambda b: spectral_radius(o, b, 1.0) - 1, 1e-3, 5.0, xtol=1e-12)
            for o in (1, 2, 3)}
    check("order 1 stability edge B_L T = 1/2", abs(edge[1] - 0.5) < 1e-9)
    x2 = edge[2] / analog_bl(2, 1.0)
    check("order 2 stability edge w0 T = sqrt6 - sqrt2",
          abs(x2 - (np.sqrt(6) - np.sqrt(2))) < 1e-9, f"{x2:.6f}")
    for f, grows in ((0.98, False), (1.02, True)):
        h = impulse_response(3, f * edge[3], 1.0, n=4000)
        check(f"order 3 impulse response {'grows' if grows else 'decays'} at"
              f" {f:.2f} x edge (B_L T = {edge[3]:.4f})",
              (np.abs(h[-200:]).max() > np.abs(h[:200]).max()) == grows)


def test_geometry():
    ps = OverheadPass()
    t = np.linspace(-200, 200, 41)
    h = 1e-3
    rho = lambda tt: ps.range_derivs(tt)[0]
    _, r1, r2, r3 = ps.range_derivs(t)
    fd1 = (rho(t + h) - rho(t - h)) / (2 * h)
    fd2 = (ps.range_derivs(t + h)[1] - ps.range_derivs(t - h)[1]) / (2 * h)
    g = 0.05
    fd3 = (ps.range_derivs(t + g)[2] - ps.range_derivs(t - g)[2]) / (2 * g)
    for name, fd, an in (("rate", fd1, r1), ("accel", fd2, r2), ("jerk", fd3, r3)):
        rel = np.max(np.abs(fd - an)) / np.max(np.abs(an))
        check(f"range {name} vs central difference", rel < 1e-6, f"rel {rel:.1e}")
    # Doppler rate at TCA by hand: rho'' = R r n^2 / (r - R)
    hand = RE * ps.r * ps.n**2 / (ps.r - RE) * ps.fc / C_LIGHT
    p1, p2, _ = ps.max_derivs()
    check("TCA Doppler rate = fc/c R r n^2/(r-R)", abs(p2 / (2 * np.pi) - hand) < 1e-6,
          f"{hand:.2f} Hz/s")
    check("max Doppler < fc v/c", p1 / (2 * np.pi) < ps.fc * np.sqrt(3.986004418e14 / ps.r) / C_LIGHT)
    check("pass ends at 10 deg elevation", abs(ps.elevation_deg(ps.t_half) - 10.0) < 1e-9)


def test_steady_state_errors():
    T = 1e-3
    t = np.arange(0, 20.0, T)
    # order 2 on a constant Doppler rate, order 3 on a constant jerk
    for order, deriv, bl in ((2, 2, 20.0), (3, 3, 5.0)):
        x = 50.0 if deriv == 2 else 5.0
        phi = x * t**deriv / (6.0 if deriv == 3 else 2.0)
        err, _ = track(order, bl, T, t, phi, np.inf)
        w0 = w0_for_bl(order, bl)
        ref = x / w0**order
        check(f"order {order} steady error = phi^({deriv}) / w0^{order}",
              abs(err[-1, 0] / ref - 1) < 1e-6, f"{err[-1, 0]:.6f} vs {ref:.6f} rad")
    t = np.arange(0, 40.0, T)
    err, _ = track(3, 5.0, T, t, 2.0 * t**2, np.inf)
    check("order 3 has zero error on a Doppler ramp",
          np.abs(err).max() > 0.01 and abs(err[-1, 0]) < 1e-9,
          f"peak {np.abs(err).max():.3f} rad -> {abs(err[-1, 0]):.1e}")


def test_pass_dynamic_error():
    ps = OverheadPass()
    T = 1e-3
    _, _, p3 = ps.max_derivs()
    t = np.arange(-ps.t_half, ps.t_half, T)
    phi, a, b, _ = ps.phase(t)
    for bl in (5.0, 8.0):
        err, _ = track(3, bl, T, t, phi, np.inf, init=(a[0], b[0]))
        qs = p3 / w0_for_bl(3, bl) ** 3
        check(f"pass max error vs quasi-static jerk bound, B_L {bl}",
              abs(np.abs(err).max() / qs - 1) < 0.02,
              f"{np.degrees(np.abs(err).max()):.3f} vs {np.degrees(qs):.3f} deg")


def test_linear_jitter():
    # High C/N0, no dynamics: phase variance -> B_L(actual) / (C/N0)
    T, cn0 = 1e-3, 45.0
    t = np.arange(0, 30.0, T)
    for order, bl in ((2, 20.0), (3, 10.0)):
        err, slips = track(order, bl, T, t, np.zeros_like(t), cn0, npaths=200, seed=11)
        var = np.var(err[5000:])
        ref = digital_bl(order, bl, T) / 10 ** (cn0 / 10)
        check(f"thermal jitter order {order}", slips.sum() == 0 and abs(var / ref - 1) < 0.03,
              f"var/ref = {var / ref:.4f}")


def test_first_order_slips():
    for rho in (2.0, 5.0, 10.0):
        q, c = mean_slip_time_quad(rho, 1.0), mean_slip_time_closed(rho, 1.0)
        check(f"slip time closed form vs Fokker-Planck quadrature rho={rho}",
              abs(q / c - 1) < 1e-8, f"B_L T_slip = {c:.4g}")
    for rho in (1.5, 2.0):
        ts = simulate_first_order_slips(rho, 10.0, npaths=2000, seed=3)
        m, ref = np.nanmean(ts), mean_slip_time_closed(rho, 10.0)
        se = np.nanstd(ts) / np.sqrt(np.count_nonzero(np.isfinite(ts)))
        check(f"Monte Carlo mean slip time rho={rho}",
              np.all(np.isfinite(ts)) and abs(m - ref) < 3.5 * se + 0.02 * ref,
              f"{m:.3f} s vs {ref:.3f} s (se {se:.3f})")
    rho = 2.0
    c = stationary_cos_mean(rho, 10.0)
    check("Tikhonov E[cos phi] = I1/I0", abs(c - i1(rho) / i0(rho)) < 0.01,
          f"{c:.4f} vs {i1(rho) / i0(rho):.4f}")


def test_optimum():
    dyn = 63.93
    for order in (2, 3):
        b = optimal_bl(order, 25.0, dyn)
        f0 = budget_deg(order, b, 25.0, dyn)
        check(f"analytic optimum is a minimum, order {order}",
              budget_deg(order, 1.01 * b, 25.0, dyn) > f0
              and budget_deg(order, 0.99 * b, 25.0, dyn) > f0)
        bd = optimal_bl(order, 25.0, dyn, T=1e-6)
        check(f"discrete optimum -> analytic as T -> 0, order {order}",
              abs(bd / b - 1) < 1e-3, f"{bd:.4f} vs {b:.4f} Hz")


if __name__ == "__main__":
    tests = [test_analog_bandwidth, test_digital_bandwidth, test_geometry,
             test_steady_state_errors, test_pass_dynamic_error, test_linear_jitter,
             test_optimum, test_first_order_slips]
    for f in tests:
        print(f.__name__)
        f()
    print("all checks passed")
