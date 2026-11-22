# Astron Built this AP
# Project 056 | Analog habitat systems / ECLSS | Multi-zone CO2 ventilation network for a sealed habitat
# Author: Aditya (aditya-astron)
# Date: 2026-11-22 IST
#
# Independent recomputation of the Rust study with numpy/scipy: the zone
# equations are rebuilt here from the parameters, integrated with an adaptive
# DOP853 solver, and compared against `cargo run --release -- --check`.

import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigvalsh, expm
from scipy.optimize import brentq, minimize_scalar

ROOT = Path(__file__).resolve().parents[1]

R, T, M, MMHG = 8.314462618, 295.15, 0.04401, 133.322368
LIMIT = 3.0 * MMHG * M / (R * T)  # kg/m^3

V = np.array([120.0, 30.0, 80.0, 60.0, 36.0])  # Hab, Node, Lab, Greenhouse, Quarters
ETA_QS = 0.85 * 60.0
AWAKE = 1.04 / (16 + 8 * 0.7)
SLEEP = 0.7 * AWAKE
S_DAY = np.array([2 * AWAKE, 0, AWAKE, AWAKE - 0.02, 0])
S_SLEEP = np.array([0, 0, 0, 0, 4 * SLEEP])


def exchange(q_quarters, eta_qs=ETA_QS):
    L = np.zeros((5, 5))
    for i, q in ((0, 150.0), (2, 100.0), (3, 80.0), (4, q_quarters)):
        L[i, 1] += q
        L[1, i] += q
        L[i, i] -= q
        L[1, 1] -= q
    L[0, 0] -= eta_qs
    return L


def rhs(L, S):
    A = L / V[:, None]
    return lambda t, c: A @ c + S / V


def integrate(L, S, c0, t, **kw):
    sol = solve_ivp(rhs(L, S), (0, t), c0, method="DOP853", rtol=1e-12, atol=1e-16, dense_output=True, **kw)
    return sol


def periodic_start(L):
    # The one-day map is affine, c -> P c + b. Recover P and b column by
    # column from ODE solves (not from the exponential formula used in Rust)
    # and solve for its fixed point.
    def day(c):
        c = integrate(L, S_DAY, c, 16.0).y[:, -1]
        return integrate(L, S_SLEEP, c, 8.0).y[:, -1]
    b = day(np.zeros(5))
    scale = 0.005
    P = np.column_stack([(day(scale * e) - b) / scale for e in np.eye(5)])
    c = np.linalg.solve(np.eye(5) - P, b)
    assert np.allclose(day(c), c, rtol=1e-10, atol=0)
    return c


def quarters_peak_mmhg(q):
    L = exchange(q)
    c0 = periodic_start(L)
    c22 = integrate(L, S_DAY, c0, 16.0).y[:, -1]
    best = 0.0
    for S, start, h in ((S_DAY, c0, 16.0), (S_SLEEP, c22, 8.0)):
        sol = integrate(L, S, start, h)
        ts = np.linspace(0, h, 961)
        ys = sol.sol(ts)[4]
        best = max(best, ys.max())
        i = int(np.argmax(ys))
        if 0 < i < len(ts) - 1:
            r = minimize_scalar(lambda t: -sol.sol(t)[4], bounds=(ts[i - 1], ts[i + 1]), method="bounded",
                                options={"xatol": 1e-10})
            best = max(best, -r.fun)
    return best * R * T / M / MMHG


def main():
    out = subprocess.run(["cargo", "run", "--release", "-q", "--", "--check"], cwd=ROOT,
                         capture_output=True, text=True, check=True).stdout
    rust = {k: float(v) for k, v in (line.split("=") for line in out.split())}

    ref = {"awake_rate": AWAKE}
    L60 = exchange(60.0)
    Vh = np.sqrt(V)
    rates = np.sort(-eigvalsh(L60 / np.outer(Vh, Vh)))
    for i, r in enumerate(rates):
        ref[f"rate{i}"] = r
    A = L60 / V[:, None]
    for name, S in (("steady_day", S_DAY), ("steady_sleep", S_SLEEP)):
        for k, v in enumerate(np.linalg.solve(A, -S / V)):
            ref[f"{name}{k}"] = v
    c0 = periodic_start(L60)
    c22 = integrate(L60, S_DAY, c0, 16.0).y[:, -1]
    for k in range(5):
        ref[f"periodic0_{k}"] = c0[k]
        ref[f"periodic2200_{k}"] = c22[k]

    # Quarters fan failure at 22:00.
    sol = integrate(exchange(3.0), S_SLEEP, c22, 8.0)
    ref["t_fan"] = brentq(lambda t: sol.sol(t)[4] - LIMIT, 0.0, 8.0, xtol=1e-13)
    # Scrubber outage at 06:00.
    sol = integrate(exchange(60.0, eta_qs=0.0), S_DAY, c0, 16.0)
    ref["t_scrub_hab"] = brentq(lambda t: sol.sol(t)[0] - LIMIT, 0.0, 16.0, xtol=1e-13)
    ref["q_min_periodic"] = brentq(lambda q: quarters_peak_mmhg(q) - 3.0, 20.0, 60.0, xtol=1e-6)

    # Also confirm the expm-based map agrees with the ODE solver directly.
    G = np.zeros((6, 6))
    G[:5, :5] = A * 16.0
    G[:5, 5] = S_DAY / V * 16.0
    via_expm = (expm(G) @ np.append(c0, 1.0))[:5]
    assert np.allclose(via_expm, c22, rtol=1e-9, atol=0)

    worst = 0.0
    tol = {"q_min_periodic": 1e-6}
    for k, v in ref.items():
        err = abs(rust[k] - v) / abs(v)
        limit = tol.get(k, 1e-8)
        flag = "ok" if err < limit else "FAIL"
        worst = max(worst, err / limit)
        print(f"{k:<18} rust {rust[k]: .10e}  scipy {v: .10e}  rel {err:.1e}  {flag}")
    if worst >= 1.0:
        sys.exit("cross-check failed")
    print(f"all {len(ref)} quantities agree")


if __name__ == "__main__":
    main()
