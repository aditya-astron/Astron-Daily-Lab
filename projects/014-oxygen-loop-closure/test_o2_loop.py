# Astron Built this AP
# Project 014 | ECLSS | Oxygen loop closure: PEM electrolysis and Sabatier CO2 reduction
# Author: Aditya (aditya-astron)
# Date: 2026-11-04 IST
"""Checks against JANAF/NIST reference data, hand calculations and the
identities the equilibrium and loop balance must satisfy.

Run: python3 test_o2_loop.py   (or pytest)
"""
import math

from thermo import (FARADAY, R, SABATIER, ELECTROLYSIS, cp, entropy,
                    reaction, ln_K)
from o2_loop import (Crew, PEMCell, loop_balance, sabatier_extent,
                     M_O2, M_H2O, DAY)


def close(a, b, tol):
    return abs(a - b) <= tol


def test_shomate_against_janaf():
    # (species, T, Cp J/mol/K, S J/mol/K) from JANAF tables
    ref = [("H2", 298.15, 28.836, 130.680), ("H2", 1000, 30.205, 166.216),
           ("O2", 298.15, 29.376, 205.147), ("CO2", 298.15, 37.135, 213.785),
           ("CO2", 1000, 54.308, 269.299), ("CH4", 298.15, 35.639, 186.251),
           ("CH4", 1000, 71.795, 247.549), ("H2O", 298.15, 33.590, 188.834),
           ("H2O", 1000, 41.268, 232.738), ("H2O(l)", 298.15, 75.351, 69.950)]
    for sp, T, c, s in ref:
        assert close(cp(sp, T), c, 0.06), (sp, T, cp(sp, T))
        assert close(entropy(sp, T), s, 0.06), (sp, T, entropy(sp, T))


def test_reaction_standard_values():
    dH, _, dG = reaction(SABATIER, 298.15)
    assert close(dH, -165.0e3, 300) and close(dG, -113.5e3, 300)
    dH, _, dG = reaction(ELECTROLYSIS, 298.15)
    assert close(dG / (2 * FARADAY), 1.229, 0.001)       # reversible voltage
    assert close(dH / (2 * FARADAY), 1.481, 0.001)       # thermoneutral voltage
    assert close(PEMCell(T=353.15).e_rev(), 1.184, 0.003)


def test_van_t_hoff_consistency():
    # d lnK / dT = dH / (R T^2) must hold if H and S come from one Cp fit.
    for T in (450.0, 600.0, 800.0):
        h = 1e-3
        num = (ln_K(SABATIER, T + h) - ln_K(SABATIER, T - h)) / (2 * h)
        ana = reaction(SABATIER, T)[0] / (R * T * T)
        assert close(num, ana, 1e-6 * abs(ana)), (T, num, ana)


def test_kirchhoff_hand_calc():
    # Independent route: 298 K formation data and a constant dCp = -49.6 J/mol/K.
    dH0, dS0, dCp, T0, T = -165.0e3, -172.6, -49.6, 298.15, 573.15
    dH = dH0 + dCp * (T - T0)
    dS = dS0 + dCp * math.log(T / T0)
    lnK_hand = -(dH - T * dS) / (R * T)
    assert close(ln_K(SABATIER, T), lnK_hand, 0.15), (ln_K(SABATIER, T), lnK_hand)


def test_equilibrium_solution():
    a, b, T, P = 1.0, 4.0, 623.15, 1.0e5
    xi = sabatier_extent(a, b, T, P)
    n = [a - xi, b - 4 * xi, xi, 2 * xi]            # CO2 H2 CH4 H2O
    tot = sum(n)
    y = [k / tot for k in n]
    lnQ = math.log(y[2]) + 2 * math.log(y[3]) - math.log(y[0]) - 4 * math.log(y[1])
    assert close(lnQ, ln_K(SABATIER, T), 1e-9)
    # element balances
    assert close(n[0] + n[2], a, 1e-12)                             # C
    assert close(2 * n[1] + 4 * n[2] + 2 * n[3], 2 * b, 1e-12)      # H
    assert close(2 * n[0] + n[3], 2 * a, 1e-12)                     # O


def test_equilibrium_trends_and_literature():
    x = [sabatier_extent(1, 4, Tc + 273.15) for Tc in range(200, 701, 50)]
    assert all(x1 > x2 for x1, x2 in zip(x, x[1:]))                 # exothermic
    assert sabatier_extent(1, 4, 773.15, 10e5) > sabatier_extent(1, 4, 773.15)   # 5 mol -> 3 mol
    # Published thermodynamic analyses, 1 atm, H2:CO2 = 4 (Gao et al. 2012)
    assert 0.93 < sabatier_extent(1, 4, 573.15) < 0.97
    assert 0.60 < sabatier_extent(1, 4, 773.15) < 0.75


def test_loop_hand_calc_and_ceiling():
    crew = Crew()
    r = loop_balance(crew, 573.15)
    n_o2 = 3.36 / M_O2
    assert close(r["water_electrolysed_kg"], 2 * n_o2 * M_H2O, 1e-9)
    assert close(r["water_electrolysed_kg"], 3.783, 0.001)
    assert close(r["stack_current_A"] * 12, 4 * FARADAY * n_o2 / DAY, 1e-9)
    assert close(r["stack_current_A"], 39.09, 0.01)
    assert close(r["water_recovered_kg"] + r["makeup_water_kg"], r["water_electrolysed_kg"], 1e-12)
    # H2-limited: closure below 1/2 at every temperature and approaching it when cold
    for Tc in (150, 250, 350, 450, 600):
        assert loop_balance(crew, Tc + 273.15)["closure"] < 0.5
    assert loop_balance(crew, 423.15)["closure"] > 0.499
    # stack energy: electrical in = thermoneutral (chemical) + waste heat
    cell = PEMCell()
    I_tot = r["stack_current_A"] * 12
    assert close(r["stack_power_W"], cell.v_tn() * I_tot + r["stack_heat_W"], 1e-9)


def test_cell_polarisation():
    cell = PEMCell()
    v = [cell.voltage(i) for i in (0.1, 0.5, 1.0, 1.5, 2.0)]
    assert all(b > a for a, b in zip(v, v[1:]))
    assert v[0] > cell.e_rev()
    assert 1.75 < cell.voltage(1.0) < 1.95       # typical PEM stack at 80 C, 1 A/cm2


if __name__ == "__main__":
    tests = [f for k, f in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS {t.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"FAIL {t.__name__}: {e}")
    print(f"{len(tests) - failed}/{len(tests)} passed")
    raise SystemExit(failed != 0)
