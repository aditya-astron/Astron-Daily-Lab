# Astron Built this AP
# Project 061 | Radiation shielding | Sector-based TID analysis and spot-shield sizing for an avionics box
# Author: Aditya (aditya-astron)
# Date: 2026-11-23 IST
#
# Independent cross-check of the Rust study. The closed-form references use
# scipy (expn for E2, quad and dblquad), and the bus model is re-traced
# with a vectorised numpy slab/sphere intersection written from scratch.
# The numbers are compared with the values printed by `cargo run --release`.

import re
import subprocess
from pathlib import Path

import numpy as np
from scipy.integrate import dblquad, quad
from scipy.special import expn

ROOT = Path(__file__).resolve().parent.parent
TERMS = [(4000.0, 0.10), (300.0, 0.40), (3.0, 5.0)]
AL = 2.70


def D(t):
    t = np.asarray(t, dtype=float)
    return sum(a * np.exp(-t / l) for a, l in TERMS)


def fib(n):
    i = np.arange(n, dtype=float)
    z = 1.0 - (2.0 * i + 1.0) / n
    r = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    phi = np.pi * (3.0 - np.sqrt(5.0)) * i
    return np.column_stack([r * np.cos(phi), r * np.sin(phi), z])


def chord_box(p, u, lo, hi):
    with np.errstate(divide="ignore", invalid="ignore"):
        ta = (np.asarray(lo) - p) / u
        tb = (np.asarray(hi) - p) / u
    t0 = np.where(np.isnan(ta), -np.inf, np.minimum(ta, tb))
    t1 = np.where(np.isnan(tb), np.inf, np.maximum(ta, tb))
    s0 = np.maximum(t0.max(axis=1), 0.0)
    s1 = t1.min(axis=1)
    return np.maximum(s1 - s0, 0.0)


def chord_sphere(p, u, c, r):
    d = p - np.asarray(c)
    b = u @ d
    disc = b * b - (d @ d - r * r)
    q = np.sqrt(np.maximum(disc, 0.0))
    s0 = np.maximum(-b - q, 0.0)
    s1 = -b + q
    return np.where(disc > 0, np.maximum(s1 - s0, 0.0), 0.0)


def baseline_bodies():
    """Same geometry as model.rs, listed directly (not imported)."""
    b = []

    def shell(lo, hi, w, rho):
        b.append(("box", lo, hi, rho))
        b.append(("box", [x + w for x in lo], [x - w for x in hi], -rho))

    shell([-80.0] * 3, [80.0] * 3, 1.0, 0.40)
    b.append(("sph", [0, 0, 0], 35.0, 4.43))
    b.append(("sph", [0, 0, 0], 34.85, -4.43))
    b.append(("sph", [0, 0, 0], 34.85, 0.50))
    b.append(("box", [20, 30, -20], [45, 60, 20], 1.20))
    b.append(("box", [40, -60, -15], [70, -30, 15], 0.60))
    shell([49.0, -10.0, -7.5], [79.0, 10.0, 7.5], 0.20, AL)
    for x in [52.0, 57.0, 62.0, 67.0, 72.0]:
        b.append(("box", [x, -9.0, -6.5], [x + 0.16, 9.0, 6.5], 1.85))
    return b


def areal(bodies, p, u):
    t = np.zeros(len(u))
    for kind, a, c, rho in bodies:
        t += rho * (chord_box(p, u, a, c) if kind == "box" else chord_sphere(p, u, a, c))
    return np.maximum(t, 0.0)


def rust_output():
    out = subprocess.run(
        ["cargo", "run", "--release", "--quiet"], cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout
    checks = {}
    for key, pat in [
        ("sphere_centre", r"sphere shell, centre\s+([\d.]+)"),
        ("sphere_offset", r"sphere shell, d = 0.6 R\s+([\d.]+)"),
        ("slab", r"two infinite plates\s+([\d.]+)"),
        ("cube", r"cubic shell, centre\s+([\d.]+)"),
    ]:
        checks[key] = float(re.search(pat, out).group(1))
    parts = {m.group(1): float(m.group(2)) for m in re.finditer(r"   (U\d) \S+ +\([^)]*\)\s+([\d.]+)\s+[\d.]+ mm", out)}
    return checks, parts


def main():
    checks, parts = rust_output()
    u = fib(200_000)

    ref = {
        "sphere_centre": float(D(0.3 * AL)),
        "sphere_offset": 0.5
        * quad(
            lambda mu: D(AL * (np.sqrt(20.0**2 - 144 + 144 * mu * mu) - np.sqrt(19.7**2 - 144 + 144 * mu * mu))),
            -1,
            1,
            epsabs=1e-13,
            epsrel=1e-13,
        )[0],
        "slab": sum(a * expn(2, 0.54 / l) for a, l in TERMS),
        "cube": 24.0
        / (4 * np.pi)
        * dblquad(
            lambda y, x: D(0.54 * np.sqrt(1 + x * x + y * y)) / (1 + x * x + y * y) ** 1.5,
            0,
            1,
            0,
            1,
            epsabs=1e-13,
            epsrel=1e-12,
        )[0],
    }
    print("closed-form references (Rust value printed to 6 decimals)")
    for k, v in ref.items():
        r = checks[k] / v - 1.0
        print(f"  {k:14s} rust {checks[k]:.6f}  scipy {v:.6f}  rel {r:+.1e}")
        assert abs(r) < 1e-6, k

    bodies = baseline_bodies()
    pts = {
        "U1": [72.26, 6.0, 4.0],
        "U2": [72.26, 0.0, 0.0],
        "U3": [62.26, 0.0, 0.0],
        "U4": [52.26, -6.0, -3.0],
    }
    print("baseline part doses, independent numpy ray trace (N = 200k)")
    for name, p in pts.items():
        d = float(D(areal(bodies, np.asarray(p), u)).mean())
        r = parts[name] / d - 1.0
        print(f"  {name}  rust {parts[name]:.2f}  numpy {d:.4f}  rel {r:+.1e}")
        assert abs(r) < 1e-3, name  # Rust prints 2 decimals
    print("all cross-checks agree")


if __name__ == "__main__":
    main()
