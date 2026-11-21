# Astron Built this AP
# Project 054 | Attitude control | Four-wheel pyramid momentum envelope and zero-crossing avoidance
# Author: Aditya (aditya-astron)
# Date: 2026-11-21 IST
"""Envelope trade over pyramid elevation and a three-orbit wheel-speed run."""
import numpy as np

from wheel_pyramid import (NULL, best_beta, fibonacci_sphere, pinv_radius,
                           pyramid, radius_facets, radius_lp, simulate,
                           sphere_radius, zero_crossings)

# Wheel: 4 N m s at 6000 rpm; stay above 600 rpm to keep clear of the
# stiction / lubricant-starvation region and the speed-sign ambiguity of the
# tachometer.
HMAX = 4.0
RPM_MAX = 6000.0
HMIN = HMAX * 600.0 / RPM_MAX
TAU_NULL = 0.02          # N m per wheel reserved for redistribution
K_NULL = 0.01            # 1/s, about 100 s time constant
PERIOD = 5700.0          # s, ~500 km LEO
N0 = 2 * np.pi / PERIOD


def body_momentum(t):
    """Stored momentum for an inertially held LEO bus: orbit-rate and twice
    orbit-rate gravity-gradient terms plus a fixed offset (N m s)."""
    return np.array([1.6 * np.sin(N0 * t),
                     1.2 * np.sin(2 * N0 * t + 0.4),
                     0.5 + 1.0 * np.cos(N0 * t)])


def rpm(h):
    return h / HMAX * RPM_MAX


def main():
    ex, ez = np.array([1.0, 0, 0]), np.array([0, 0, 1.0])
    print("Envelope vs pyramid elevation (units of hmax)")
    print(" beta[deg]  r_x(LP)  2r2*cos  r_z(LP)  4*sin   sphere")
    for bd in (20.0, 30.0, np.rad2deg(np.arctan(1 / np.sqrt(2))), 45.0, 60.0):
        b = np.deg2rad(bd)
        A = pyramid(b)
        print(f" {bd:8.3f}  {radius_lp(A, ex):7.4f}  {2*np.sqrt(2)*np.cos(b):7.4f}"
              f"  {radius_lp(A, ez):7.4f}  {4*np.sin(b):6.4f}  {sphere_radius(A):6.4f}")

    b, r = best_beta()
    print(f"\nbest elevation {np.rad2deg(b):.4f} deg (atan(1/sqrt2) = "
          f"{np.rad2deg(np.arctan(1/np.sqrt(2))):.4f}), sphere {r:.6f} hmax "
          f"(4/sqrt6 = {4/np.sqrt(6):.6f})")

    A = pyramid(b)
    D = fibonacci_sphere(5000)
    env = np.array([radius_facets(A, d) for d in D])
    pin = np.array([pinv_radius(A, d) for d in D])
    i = np.argmin(pin / env)
    a1 = A[:, 0]
    print(f"along wheel 1 axis: envelope {radius_lp(A, a1):.4f}, min-norm "
          f"{pinv_radius(A, a1):.4f} hmax (closed forms 2 and 4/3)")
    print(f"envelope over 5000 directions: {env.min():.4f} .. {env.max():.4f} hmax")
    print(f"min-norm law alone reaches {(pin/env).min():.4f} of the envelope at "
          f"d = {np.round(D[i], 3)}")
    scale = (HMAX - HMIN) / (2 * HMAX)
    print(f"zero-crossing-free envelope = {scale:.3f} x full envelope; "
          f"guaranteed {scale * r * HMAX:.3f} N m s in any direction")

    print(f"\nThree-orbit run, wheels {HMAX} N m s at {RPM_MAX:.0f} rpm, "
          f"floor {rpm(HMIN):.0f} rpm, null torque {TAU_NULL} N m")
    for label, use in (("min-norm only", False), ("null-motion bias", True)):
        t, W, H, feas = simulate(A, body_momentum, 3 * PERIOD, 1.0, HMIN, HMAX,
                                 TAU_NULL, K_NULL, use_null=use)
        settled = t > 600.0
        low = np.mean(np.abs(W[settled]) < HMIN, axis=0)
        err = np.max(np.abs(W @ A.T - H))
        print(f"  {label:17s} zero crossings {zero_crossings(W[settled]):3d}  "
              f"min |w| {rpm(np.abs(W[settled]).min()):6.1f} rpm  "
              f"max |w| {rpm(np.abs(W).max()):6.1f} rpm  "
              f"worst wheel below floor {100*low.max():5.1f} %  |Aw-h| {err:.1e}")
    print(f"  null vector {NULL}, |A n| = {np.linalg.norm(A @ NULL):.1e}")


if __name__ == "__main__":
    main()
