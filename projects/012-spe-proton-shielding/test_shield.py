# Astron Built this AP
# Project 012 | Radiation shielding | Solar proton event dose behind aluminium vs polyethylene shielding
# Author: Aditya (aditya-astron)
# Date: 2026-11-04 IST
"""Checks against NIST PSTAR and against identities the transport must satisfy.

Run: python3 test_shield.py
"""
import numpy as np

from stopping import E_MIN, WATER, ALUMINIUM, POLYETHYLENE
from shield import (MEV_PER_G_TO_GY, king_aug1972, king_integral, transmit,
                    spectrum_out, dose_normal, dose_isotropic)

failures = 0


def check(name, ok, detail=""):
    global failures
    failures += not ok
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {detail}")


def rel(a, b):
    return abs(a - b) / abs(b)


# NIST PSTAR, liquid water: E (MeV), S_electronic (MeV cm^2/g), CSDA range (g/cm^2)
PSTAR_WATER = [(10, 45.67, 0.1230), (100, 7.289, 7.718), (200, 4.492, 25.96)]
for E, S_ref, R_ref in PSTAR_WATER:
    S, R = float(WATER.S(E)), float(WATER.range(E))
    check(f"water S({E} MeV)", rel(S, S_ref) < 0.01, f"{S:.3f} vs {S_ref}")
    check(f"water CSDA R({E} MeV)", rel(R, R_ref) < 0.015, f"{R:.4f} vs {R_ref}")

# Range table inverts cleanly.
E = np.geomspace(2, 2000, 50)
for m in (WATER, ALUMINIUM, POLYETHYLENE):
    err = np.max(np.abs(m.energy(m.range(E)) / E - 1))
    check(f"{m.name} E(R(E)) round trip", err < 1e-9, f"max rel err {err:.1e}")

# Per gram, hydrogen-rich polyethylene out-stops aluminium (higher Z/A, lower I).
check("HDPE stops more than Al per g/cm^2",
      bool(np.all(POLYETHYLENE.S(E) > ALUMINIUM.S(E))))

# Zero thickness is the identity.
E_in, jac = transmit(E, [])
check("empty stack is identity", np.allclose(E_in, E) and np.allclose(jac, 1))

# Two layers of the same material equal one of the summed thickness.
a, _ = transmit(E, [(ALUMINIUM, 3.0), (ALUMINIUM, 4.0)])
b, _ = transmit(E, [(ALUMINIUM, 7.0)])
check("3 + 4 g/cm^2 Al == 7 g/cm^2 Al", np.allclose(a, b, rtol=1e-9))

# Proton conservation: everything leaving above E_MIN entered above E(R(E_MIN)+t).
for m, t in [(ALUMINIUM, 2.0), (POLYETHYLENE, 10.0)]:
    Eo = np.geomspace(E_MIN, 2999.0, 20000)
    n_out = np.trapezoid(spectrum_out(king_aug1972, [(m, t)], Eo) * Eo, np.log(Eo))
    E_c = float(m.energy(m.range(E_MIN) + t))
    n_ref = king_integral(E_c)
    check(f"proton count through {t} g/cm^2 {m.name}", rel(n_out, n_ref) < 2e-3,
          f"{n_out:.4e} vs {n_ref:.4e}")

# Unshielded dose has a closed form for the exponential spectrum only numerically,
# but its integrand is analytic: compare against a fine independent quadrature.
from scipy.integrate import quad
d_ref = MEV_PER_G_TO_GY * quad(lambda e: king_aug1972(e) * WATER.S(e), 10, 3000,
                               limit=400, points=[30, 100, 300])[0]
d0 = dose_normal(king_aug1972, [])
check("unshielded dose vs scipy quad", rel(d0, d_ref) < 1e-3, f"{d0:.3f} vs {d_ref:.3f} Gy")

# Isotropic exposure through nothing: half the omnidirectional free-space dose.
check("D_iso(0) = D_free / 2", rel(dose_isotropic(king_aug1972, ALUMINIUM, 0.0), 0.5 * d0) < 1e-12)

# Thin slab: isotropic result converges to the zero-thickness limit.
d_thin = dose_isotropic(king_aug1972, ALUMINIUM, 1e-4)
check("D_iso(t -> 0) continuous", rel(d_thin, 0.5 * d0) < 0.01, f"{d_thin:.3f} vs {0.5 * d0:.3f}")

# Monoenergetic beam: 150 MeV through 10 g/cm^2 water leaves at E(R(150) - 10).
E_exit = float(WATER.energy(WATER.range(150.0) - 10.0))
E_back, _ = transmit(np.array([E_exit]), [(WATER, 10.0)])
check("150 MeV through 10 g/cm^2 water", rel(E_back[0], 150.0) < 1e-9, f"exit {E_exit:.2f} MeV")

# Dose falls monotonically with thickness, and HDPE beats Al at every depth.
ts = [0.5, 1, 2, 5, 10, 20]
dal = [dose_isotropic(king_aug1972, ALUMINIUM, t) for t in ts]
dpe = [dose_isotropic(king_aug1972, POLYETHYLENE, t) for t in ts]
check("dose decreases with thickness", bool(np.all(np.diff(dal) < 0) and np.all(np.diff(dpe) < 0)))
check("HDPE dose < Al dose at equal areal density", all(p < a for p, a in zip(dpe, dal)))

# Gauss-Legendre in mu has converged.
d48 = dose_isotropic(king_aug1972, ALUMINIUM, 5.0, n_mu=48)
d96 = dose_isotropic(king_aug1972, ALUMINIUM, 5.0, n_mu=96)
check("mu quadrature converged at 5 g/cm^2", rel(d48, d96) < 1e-3, f"{d48:.4f} vs {d96:.4f}")

print(f"\n{failures} failure(s)")
raise SystemExit(failures != 0)
