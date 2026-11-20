# Astron Built this AP
# Project 050 | Rocket propulsion | Regenerative cooling of a LOX/methane thrust chamber
# Author: Aditya (aditya-astron)
# Date: 2026-11-20 IST
"""Checks against tables, closed forms and an independent 2D conduction solve."""
import math

import numpy as np
import pytest

from regen import (CombustionGas, ConstantFluid, CoolPropFluid, Jacket,
                   RegenSolver, area_ratio, bartz_gas_side, bartz_htc, bartz_sigma,
                   mach_from_area, solve_station)
from run_regen import CONTOUR, GAS, P_IN, T_IN, baseline_jacket, engine_flows
from unitcell import fin_model_flat, solve_cell


# --- isentropic flow -----------------------------------------------------------
def test_area_mach_against_naca_1135():
    # NACA 1135, gamma = 1.4: A/A* = 2.0 at M = 0.3059 and M = 2.1972
    assert mach_from_area(2.0, 1.4, False) == pytest.approx(0.30590, abs=5e-5)
    assert mach_from_area(2.0, 1.4, True) == pytest.approx(2.19719, abs=5e-5)
    assert area_ratio(3.0, 1.4) == pytest.approx(4.2346, abs=1e-4)


@pytest.mark.parametrize("M", [0.05, 0.3, 0.8, 1.2, 2.5, 4.0])
def test_area_mach_round_trip(M):
    eps = area_ratio(M, 1.15)
    assert mach_from_area(eps, 1.15, M > 1) == pytest.approx(M, rel=1e-10)


def test_recovery_temperature_limits():
    g = CombustionGas(3550.0, 1.15, 21.0, 7e6)
    assert g.recovery_temperature(0.0) == pytest.approx(3550.0)
    r = g.prandtl ** (1 / 3)
    M = 2.0
    s = 0.075 * M * M
    assert g.recovery_temperature(M) == pytest.approx(3550.0 * (1 + r * s) / (1 + s))
    assert g.recovery_temperature(M) < 3550.0


def test_cstar_hand_value():
    # sqrt(R T0 / g) / (2/(g+1))^((g+1)/(2(g-1))) with R = 8314.46/21
    R = 8314.462618 / 21.0
    expect = math.sqrt(R * 3550.0 / 1.15) / (2 / 2.15) ** (2.15 / 0.30)
    assert GAS.cstar == pytest.approx(expect, rel=1e-12)
    assert 1800.0 < GAS.cstar < 1900.0


# --- Bartz in SI vs. the original US customary form ----------------------------
def test_bartz_si_matches_us_customary():
    Dt, rc, eps, M, Twg = 0.06, 0.0282, 1.0, 1.0, 850.0
    h_si = bartz_htc(GAS, Dt, rc, eps, M, Twg)
    # Huzel & Huang form: Dt [in], mu [lb/(in s)], cp [Btu/(lb F)], pc [psia],
    # g = 32.174 ft/s^2, c* [ft/s]; result in Btu/(in^2 s F)
    Dt_in, rc_in = Dt / 0.0254, rc / 0.0254
    mu = 46.6e-10 * 21.0 ** 0.5 * (3550.0 * 1.8) ** 0.6
    cp = GAS.cp / 4186.8
    pc = GAS.pc / 6894.757
    cs = GAS.cstar / 0.3048
    sigma = bartz_sigma(Twg, GAS.T0, M, GAS.gamma)
    h_us = (0.026 / Dt_in ** 0.2) * (mu ** 0.2 * cp / GAS.prandtl ** 0.6) \
        * (pc * 32.174 / cs) ** 0.8 * (Dt_in / rc_in) ** 0.1 * eps ** -0.9 * sigma
    to_si = 1055.056 / (0.0254 ** 2 * (5.0 / 9.0))
    assert h_si == pytest.approx(h_us * to_si, rel=2e-4)


def test_bartz_sigma_unity_at_stagnation_and_adiabatic_wall():
    assert bartz_sigma(3550.0, 3550.0, 0.0, 1.15) == pytest.approx(1.0)


# --- contour -------------------------------------------------------------------
def test_contour_continuity_and_lstar():
    c = CONTOUR
    for xj in (c.xf, c.x2, c.x1, 0.0, c.x3):
        d = 1e-9
        assert c.r(xj - d) == pytest.approx(c.r(xj + d), abs=1e-7)
        assert c.drdx(xj - d) == pytest.approx(c.drdx(xj + d), abs=1e-4)
    assert c.r(0.0) == pytest.approx(c.Rt)
    assert c.r(c.x_exit) == pytest.approx(c.Re)
    vol = c._volume(c.x_inj, 0.0, n=20000)
    assert vol / (math.pi * c.Rt ** 2) == pytest.approx(c.lstar, rel=1e-5)


# --- march vs. closed form -----------------------------------------------------
class Cylinder:
    def __init__(self, R, L):
        self.R, self.L = R, L
        self.x_start, self.x_end = 0.0, L

    def r(self, x):
        return self.R

    def drdx(self, x):
        return 0.0


def constant_case(n):
    geo = Cylinder(0.04, 0.3)
    jk = Jacket(n_channels=60, height=2.5e-3, land=1.0e-3, t_wall=1.0e-3, k_wall=300.0)
    fluid = ConstantFluid(cp=4000.0, rho=400.0, mu=1.0e-4, k=0.15)
    hg, Taw = 15e3, 3200.0
    mdot = 1.5
    sol = RegenSolver(geo, jk, fluid, mdot, lambda x, Twg: (hg, Taw, 0.2)).run(150.0, 10e6, n=n)
    # per-channel series resistance per unit length
    w, b, h, area, Dh = jk.channel(geo.R)
    st = fluid.state(fluid.cp * 150.0, 10e6)
    Re = (mdot / 60 / area) * Dh / st.mu
    Pr = st.mu * st.cp / st.k
    x = -1.8 * math.log10((jk.roughness / Dh / 3.7) ** 1.11 + 6.9 / Re)
    f = 1 / x ** 2
    Nu = (f / 8) * (Re - 1000) * Pr / (1 + 12.7 * math.sqrt(f / 8) * (Pr ** (2 / 3) - 1))
    hc = Nu * st.k / Dh
    eta = math.tanh(math.sqrt(2 * hc / (300.0 * b)) * h) / (math.sqrt(2 * hc / (300.0 * b)) * h)
    P = 2 * math.pi * geo.R / 60
    Rp = 1 / (hg * P) + 1.0e-3 / (300.0 * P) + 1 / (hc * (w + 2 * eta * h))
    Ttot = 60 / Rp  # W/(m K) for the whole jacket
    T_out = Taw - (Taw - 150.0) * math.exp(-Ttot * geo.L / (mdot * fluid.cp))
    G = mdot / 60 / area
    dp = f * geo.L / Dh * G * G / (2 * st.rho)
    return sol, T_out, dp


def test_constant_property_march_matches_exponential():
    sol, T_out, dp = constant_case(200)
    assert sol.Tb[-1] == pytest.approx(T_out, abs=1e-3)
    assert sol.dp == pytest.approx(dp, rel=1e-9)
    assert sol.Q_total == pytest.approx(1.5 * 4000.0 * (sol.Tb[-1] - 150.0), rel=1e-12)


def test_heun_march_is_second_order():
    errs = []
    for n in (5, 10, 20):
        sol, T_out, _ = constant_case(n)
        errs.append(abs(sol.Tb[-1] - T_out))
    assert errs[0] / errs[1] == pytest.approx(4.0, rel=0.05)
    assert errs[1] / errs[2] == pytest.approx(4.0, rel=0.05)


# --- real chamber --------------------------------------------------------------
@pytest.fixture(scope="module")
def chamber():
    fluid = CoolPropFluid("Methane")
    _, mdot_f = engine_flows()
    gs = bartz_gas_side(GAS, CONTOUR)
    jk = baseline_jacket()
    return fluid, gs, jk, mdot_f, RegenSolver(CONTOUR, jk, fluid, mdot_f, gs)


def test_energy_balance_against_recorded_flux(chamber):
    fluid, gs, jk, mdot_f, solver = chamber
    s = solver.run(T_IN, P_IN, n=400)
    ds = np.sqrt(1 + np.array([CONTOUR.drdx(x) for x in s.x]) ** 2)
    integrand = s.q * 2 * np.pi * s.r * ds
    Q_quad = abs(np.trapezoid(integrand, s.x))
    assert s.Q_total == pytest.approx(Q_quad, rel=2e-3)
    assert s.Q_total == pytest.approx(mdot_f * (s.h[-1] - s.h[0]), rel=1e-12)
    # heat flux peaks just upstream of the throat, as Bartz predicts
    assert abs(s.x[np.argmax(s.q)]) < 0.005


def test_grid_convergence(chamber):
    *_, solver = chamber
    T = [solver.run(T_IN, P_IN, n=n).Tb[-1] for n in (100, 200, 400)]
    assert abs(T[2] - T[1]) < 0.4 * abs(T[1] - T[0])
    assert abs(T[2] - T[1]) < 0.5


def test_station_balance_closes(chamber):
    fluid, gs, jk, mdot_f, _ = chamber
    st = fluid.state(fluid.h_from_Tp(290.0, 11.8e6), 11.8e6)
    muw = lambda T: fluid.viscosity(T, 11.8e6)
    hg_fun = lambda Twg: gs(0.0, Twg)[0]
    Taw = gs(0.0, 1000.0)[1]
    res = solve_station(jk, CONTOUR.Rt, hg_fun, Taw, st, mdot_f / jk.n_channels, muw)
    w, b, h, _, _ = jk.channel(CONTOUR.Rt)
    q_in = res.q * 2 * math.pi * CONTOUR.Rt / jk.n_channels
    q_out = res.hc * (res.Twc - st.T) * (w + 2 * res.eta * h)
    assert q_in == pytest.approx(q_out, rel=1e-8)
    assert res.Twg - res.Twc == pytest.approx(res.q * jk.t_wall / jk.k_wall, rel=1e-12)
    assert res.hg == pytest.approx(hg_fun(res.Twg))


# --- 2D unit cell ----------------------------------------------------------------
def test_unit_cell_slab_limit_is_exact():
    # No land: the cell is a slab cooled across its full width, T linear in y.
    w, tw, k, hg, Taw, hc, Tb = 2e-3, 1e-3, 300.0, 2e4, 3000.0, 5e4, 300.0
    Q, Ts, _, _ = solve_cell(w, 0.0, tw, 0.0, k, hg, Taw, hc, Tb, n=8)
    q = (Taw - Tb) / (1 / hg + tw / k + 1 / hc)
    assert Q == pytest.approx(q * w / 2, rel=1e-10)
    assert np.allclose(Ts, Taw - q / hg, atol=1e-8)


def test_unit_cell_mesh_convergence_and_fin_model():
    cell = (0.81e-3, 0.8e-3, 0.7e-3, 3.0e-3, 320.0, 2.4e4, 3400.0, 6.7e4, 290.0)
    Q = [solve_cell(*cell, n=n)[0] for n in (10, 20, 40)]
    assert abs(Q[2] - Q[1]) < 0.5 * abs(Q[1] - Q[0])
    Q1, T1 = fin_model_flat(*cell)
    assert Q[2] == pytest.approx(Q1, rel=0.02)
    # the 1D model averages the hot face; the 2D peak sits above it
    assert solve_cell(*cell, n=40)[1].max() > T1
