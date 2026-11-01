# Astron Built this AP
# Project 004 | Spacecraft power systems | Solar array MPPT through a LEO eclipse cycle
# Author: Aditya (aditya-astron)
# Date: 2026-11-01 IST
"""Checks against the datasheet, independent solvers and hand calculations."""

import math

import numpy as np
import pytest

from mppt_sim import (SIGMA, IncrementalConductance, PerturbObserve, Scenario,
                      eclipse_half_angle, irradiance_profile, orbit_period,
                      panel_temperature, run, summarize)
from pv_model import (G_AM0, Array, CellDatasheet, cell_current,
                      cell_current_newton, extract_params,
                      open_circuit_voltage, params_at, true_mpp)

DS = CellDatasheet()
REF = extract_params(DS)


def test_fit_hits_all_four_datasheet_conditions():
    i = cell_current([0.0, DS.vmp, DS.voc], REF)
    assert i[0] == pytest.approx(DS.isc, abs=1e-10)
    assert i[1] == pytest.approx(DS.imp, abs=1e-10)
    assert abs(i[2]) < 1e-10
    # dP/dV = 0 at the datasheet MPP, by central difference
    h = 1e-6
    p = lambda v: v * float(cell_current(v, REF))
    assert abs((p(DS.vmp + h) - p(DS.vmp - h)) / (2 * h)) < 1e-6


def test_fitted_resistances_are_physical():
    assert 0.0 < REF.rs < 0.2         # tens of mOhm for a 30 cm^2 cell
    assert REF.rsh > 500.0
    assert REF.iph >= DS.isc          # Rs/Rsh drop puts Iph slightly above Isc


def test_hand_calcs_efficiency_and_fill_factor():
    pmp = DS.vmp * DS.imp
    assert pmp / (DS.voc * DS.isc) == pytest.approx(0.8643, abs=1e-4)
    assert pmp / (G_AM0 * DS.area_m2) == pytest.approx(0.2939, abs=1e-4)


@pytest.mark.parametrize("t_c", [-90.0, 28.0, 80.0])
def test_lambert_w_form_matches_newton(t_c):
    p = params_at(REF, DS, 0.7 * G_AM0, t_c)
    voc = open_circuit_voltage(p)
    for v in np.linspace(0.0, voc, 25):
        assert float(cell_current(v, p)) == pytest.approx(
            cell_current_newton(v, p), abs=1e-12)


def test_golden_section_mpp_matches_dense_sweep():
    arr = Array()
    conds = [params_at(REF, DS, g, t) for g, t in
             [(G_AM0, -80.0), (G_AM0, 40.0), (0.3 * G_AM0, 0.0)]]
    v_opt, p_opt = true_mpp(arr, conds)
    for c, vo, po in zip(conds, v_opt, p_opt):
        v = np.linspace(0.0, arr.voc(c), 400001)
        p = arr.power(v, c)
        assert po == pytest.approx(p.max(), rel=1e-9)
        assert vo == pytest.approx(v[p.argmax()], abs=1e-3)


def test_temperature_coefficients():
    hot = params_at(REF, DS, G_AM0, DS.t_ref_c + 30.0)
    cold = params_at(REF, DS, G_AM0, DS.t_ref_c - 30.0)
    # Voc follows the datasheet coefficient by construction
    dvoc = (open_circuit_voltage(hot) - open_circuit_voltage(cold)) / 60.0
    assert dvoc == pytest.approx(DS.dv_oc_dt, rel=1e-9)
    disc = (float(cell_current(0.0, hot)) - float(cell_current(0.0, cold))) / 60.0
    assert disc == pytest.approx(DS.di_sc_dt, rel=0.01)
    # Pmp coefficient is an output, not an input. Triple-junction GaAs is
    # typically -0.2 to -0.3 %/K.
    one = Array(1, 1)
    _, p = true_mpp(one, [hot, cold])
    gamma = (p[0] - p[1]) / 60.0 / (DS.vmp * DS.imp)
    assert -0.0030 < gamma < -0.0018


def test_low_irradiance_scales_current_not_voltage():
    half = params_at(REF, DS, 0.5 * G_AM0, DS.t_ref_c)
    assert float(cell_current(0.0, half)) == pytest.approx(DS.isc / 2, rel=1e-3)
    # Voc drops by about a ln(2) per thermal voltage
    dv = DS.voc - open_circuit_voltage(half)
    assert dv == pytest.approx(REF.a * math.log(2.0), rel=0.02)


def test_orbit_geometry_by_hand():
    r = 6378.137e3 + 500e3
    t_hand = 2 * math.pi * math.sqrt(r ** 3 / 3.986004418e14)
    assert orbit_period(500e3) == pytest.approx(t_hand, rel=1e-12)
    assert t_hand / 60 == pytest.approx(94.62, abs=0.01)
    # beta = 0: shadow half-angle is asin(Re / r)
    assert eclipse_half_angle(500e3, 0.0) == pytest.approx(
        math.asin(6378.137e3 / r), rel=1e-12)
    # high beta: no eclipse at all
    assert eclipse_half_angle(500e3, 75.0) == 0.0


def test_irradiance_profile_eclipse_duration():
    sc = Scenario(n_orbits=1)
    t, g, period, t_ecl = irradiance_profile(sc)
    dark = (g == 0).sum() * sc.dt
    # symmetric ramps of penumbra_s/2 shave that off the umbra on each side
    assert dark == pytest.approx(t_ecl - sc.penumbra_s, abs=2 * sc.dt)
    assert g.max() == sc.g_sun


def test_panel_equilibrium_matches_radiative_balance():
    sc = Scenario()
    eta = 0.29
    g = np.full(200000, sc.g_sun)     # 5.5 h of constant sun
    temps = panel_temperature(g, sc, eff_elec=eta)
    q = (sc.absorptance - eta) * sc.g_sun + sc.q_earth_ir
    t_eq = (q / ((sc.emit_front + sc.emit_back) * SIGMA)) ** 0.25 - 273.15
    assert temps[-1] == pytest.approx(t_eq, abs=0.05)


@pytest.mark.parametrize("cls", [PerturbObserve, IncrementalConductance])
def test_trackers_converge_on_static_curve(cls):
    arr = Array()
    c = params_at(REF, DS, G_AM0, 20.0)
    v_star, p_star = true_mpp(arr, [c])
    trk = cls(30.0, step=0.1)
    v = 30.0
    powers = []
    for _ in range(600):
        i = float(arr.current(v, c))
        powers.append(v * i)
        v = trk.update(v, i)
    assert abs(v - v_star[0]) < 0.5
    assert np.mean(powers[-100:]) > 0.999 * p_star[0]


def test_full_orbit_energy_ranking():
    res = run(Scenario(n_orbits=2))
    eff = summarize(res)
    assert eff["P&O"] > 0.995
    assert eff["IncCond"] > 0.995
    # DET at a fixed 35 V cannot use the cold-panel voltage margin
    assert eff["DET"] < 0.96
    for p in res["p_out"].values():
        assert np.all(p <= res["p_mpp"] + 1e-9)
