# Astron Built this AP
# Project 066 | Electronics and circuit analysis | Hot-swap inrush control and MOSFET junction temperature
# Author: Aditya (aditya-astron)
# Date: 2026-11-26 IST
"""Checks against hand calculation, the ramp closed forms, the exact Foster
step response and an independent ODE solve. Run `pytest -q`."""
import numpy as np
import pytest
from scipy.integrate import solve_ivp

import hotswap as h

HS = h.HotSwap(p_conv=120.0)


def test_fet_model_continuous_at_saturation_edge():
    f = h.Fet()
    vgs = 5.0
    vov = vgs - f.vth
    assert f.id(vgs, vov) == pytest.approx(f.id(vgs, vov + 1e-9), rel=1e-8)
    assert f.id(vgs, vov - 1e-9) == pytest.approx(0.5 * f.k * vov ** 2, rel=1e-8)
    assert f.rds_on(15.0) == pytest.approx(1 / (20.0 * 12.0))


def test_hand_values():
    # S = 1 V/ms: Cg = 10 uA / 1000 V/s = 10 nF, I = 0.47 A, E = 0.5*470u*28^2
    assert HS.cg_for_slew(1000.0) == pytest.approx(10e-9)
    assert h.energy_capacitive(HS) == pytest.approx(0.18424, rel=1e-12)
    # E_R = 28^3 / (6 * 50 * 1000) = 0.0731733
    assert h.energy_resistive(h.with_(HS, r_load=50.0), 1000.0) == pytest.approx(0.0731733, rel=1e-6)
    # E_P = 0.12 * (28 * 0.441833 - 10) = 0.284558
    assert h.energy_converter(HS, 1000.0) == pytest.approx(0.284558, rel=1e-5)


@pytest.mark.parametrize("slew", [100.0, 1000.0, 20000.0])
def test_capacitive_energy_exact_for_any_slew(slew):
    x = h.with_(HS, p_conv=0.0, r_load=1e12)
    tr = h.simulate(x, slew)
    assert tr.vo[-1] == pytest.approx(x.vb, abs=1e-6)
    # Cgs feeds a few nC into the output that the FET does not carry
    assert tr.energy == pytest.approx(h.energy_capacitive(x), rel=2e-4)


@pytest.mark.parametrize("slew", [300.0, 3000.0])
def test_resistive_and_power_good_closed_forms(slew):
    x = h.with_(HS, p_conv=0.0, r_load=50.0)
    assert h.simulate(x, slew).energy == pytest.approx(h.energy_closed_form(x, slew), rel=5e-3)
    pg = h.with_(HS, pg_gated=True)
    assert h.simulate(pg, slew).energy == pytest.approx(h.energy_closed_form(pg, slew), rel=5e-3)


def test_uvlo_stall_explains_excess_at_slow_slew():
    for s in (100.0, 300.0):
        tr = h.simulate(HS, s)
        excess = tr.energy - h.energy_closed_form(HS, s)
        assert excess > 0
        assert h.energy_turn_on_stall(HS, s) == pytest.approx(excess, rel=0.12)


def test_output_dip_inside_uvlo_hysteresis():
    # converter assumed latched on: the dip after turn-on must stay well under
    # a typical 2 V UVLO hysteresis
    for s in (100.0, 1000.0, 10000.0):
        tr = h.simulate(HS, s)
        assert HS.v_en - tr.vo[tr.t >= tr.t_enable].min() < 1.0


@pytest.mark.parametrize("slew", [500.0, 5000.0])
def test_peak_current_closed_form(slew):
    for pg in (False, True):
        x = h.with_(HS, pg_gated=pg)
        assert h.simulate(x, slew).id.max() == pytest.approx(
            h.peak_current_closed_form(x, slew), rel=0.07)


def test_slew_and_ramp_time():
    s = 2000.0
    tr = h.simulate(h.with_(HS, p_conv=0.0), s)
    mid = (tr.vo > 5.0) & (tr.vo < 25.0)
    rate = np.polyfit(tr.t[mid], tr.vo[mid], 1)[0]
    assert rate == pytest.approx(s, rel=1e-3)
    # 14 ms of following at Ig/Cg plus 12 V of still output at Ig/(Cg+Cgs)
    assert HS.gate_clamp_time(HS.cg_for_slew(s)) == pytest.approx(0.014 + 12 * 8e-9 / 10e-6)
    assert tr.t_done == pytest.approx(HS.gate_clamp_time(HS.cg_for_slew(s)), rel=0.01)


def test_foster_step_matches_zth():
    f = HS.fet
    for tp in (3e-5, 1e-3, 0.05, 2.0):
        t = np.array([0.0, 1e-15, tp])
        _, rise = h.junction_rise(f, t, np.array([0.0, 7.0, 7.0]), t_cool=0)
        assert rise[-1] == pytest.approx(7.0 * h.zth(f, tp)[0], rel=1e-9)
    assert h.zth(f, 1e6)[0] == pytest.approx(sum(f.r_th))


def test_foster_stepper_matches_ode_on_ramp_waveform():
    tr = h.simulate(HS, 1000.0)
    t, rise = h.junction_rise(HS.fet, tr.t, tr.p, t_cool=0.05)
    f = HS.fet
    r, tau = np.array(f.r_th), np.array(f.tau_th)
    # piecewise-linear power, same as the stepper assumes
    pfun = lambda s: np.interp(s, tr.t, tr.p, right=0.0)
    sol = solve_ivp(lambda s, th: (r * pfun(s) - th) / tau, (0.0, t[-1]), np.zeros(len(r)),
                    method="Radau", rtol=1e-9, atol=1e-12, t_eval=t[::50],
                    max_step=2e-5)
    assert np.max(np.abs(sol.y.sum(axis=0) - rise[::50])) < 1e-5


def test_peak_tj_converged_in_sampling():
    a, _ = h.peak_tj(HS, 2000.0, n_out=4001)
    b, _ = h.peak_tj(HS, 2000.0, n_out=40001)
    assert a == pytest.approx(b, abs=0.05)


def test_design_trends():
    pg = h.with_(HS, pg_gated=True)
    tj_pg = [h.peak_tj(pg, s, n_out=3001)[0] for s in (200.0, 1000.0, 5000.0)]
    assert tj_pg[0] < tj_pg[1] < tj_pg[2]           # monotonic: slower is cooler
    tj_uv = [h.peak_tj(HS, s, n_out=3001)[0] for s in (200.0, 2000.0, 20000.0)]
    assert tj_uv[1] < tj_uv[0] and tj_uv[1] < tj_uv[2]   # interior optimum
    # power-good gating: ramp heating independent of converter power
    a = h.peak_tj(h.with_(pg, p_conv=60.0), 1000.0, n_out=3001)[0]
    b = h.peak_tj(h.with_(pg, p_conv=600.0), 1000.0, n_out=3001)[0]
    assert a == pytest.approx(b, abs=0.1)
