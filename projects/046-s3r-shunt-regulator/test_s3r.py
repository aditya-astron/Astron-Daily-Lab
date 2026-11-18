# Astron Built this AP
# Project 046 | Spacecraft power systems | S3R shunt regulator for a 50 V GEO bus
# Author: Aditya (aditya-astron)
# Date: 2026-11-18 IST
import math

import numpy as np
import pytest
from scipy.optimize import brentq

from s3r import (CellSpec, SectionDesign, S3R, make_cell, max_power_point,
                 measure_cycle, steady_state_theory)

SPEC = CellSpec()


# ------------------------------------------------------------------ array

@pytest.mark.parametrize("temp,eol", [(-20, False), (28, False), (70, True)])
def test_cell_pinned_to_isc_and_voc(temp, eol):
    cell = make_cell(SPEC, temp, eol)
    isc = SPEC.isc + SPEC.disc_dt * (temp - SPEC.t_ref)
    voc = SPEC.voc + SPEC.dvoc_dt * (temp - SPEC.t_ref)
    if eol:
        isc, voc = isc * SPEC.eol_isc, voc * SPEC.eol_voc
    assert cell.current(0.0) == pytest.approx(isc, rel=1e-12)
    assert abs(cell.current(voc)) < 1e-12


def test_lambert_w_matches_bracketed_root():
    cell = make_cell(SPEC, 55.0, eol=True, sun=0.9)
    for v in np.linspace(0.0, 2.6, 27):
        for rs_extra in (0.0, 0.03):
            i_ref = brentq(lambda i: cell.residual(v, i, rs_extra),
                           -5.0, cell.iph + cell.i0, xtol=1e-15, rtol=1e-15)
            assert cell.current(v, rs_extra) == pytest.approx(i_ref, abs=1e-12)


def test_large_argument_branch_continuous():
    # Far past Voc the W argument overflows exp(); the asymptotic branch
    # must agree with the root of the implicit equation.
    cell = make_cell(SPEC, 28.0, eol=False)
    v = 25.0
    i_ref = brentq(lambda i: cell.residual(v, i), -1e3, 0.0, xtol=1e-14)
    assert cell.current(v) == pytest.approx(i_ref, rel=1e-10)


def test_cell_mpp_is_realistic():
    # Commercial 3J GaAs cells of this size: Vmp ~2.4 V, FF ~0.85.
    vmp, imp = max_power_point(make_cell(SPEC, 28.0, eol=False))
    assert 2.30 < vmp < 2.45
    assert 0.83 < vmp * imp / (SPEC.voc * SPEC.isc) < 0.88


def test_irradiance_scales_short_circuit_current():
    a = make_cell(SPEC, 40.0, eol=True, sun=1.0)
    b = make_cell(SPEC, 40.0, eol=True, sun=0.8876)
    # Rs drop at short circuit makes this very nearly, not exactly, linear.
    assert b.current(0.0) / a.current(0.0) == pytest.approx(0.8876, rel=1e-6)


def test_section_is_stiff_current_source_at_bus():
    sd = SectionDesign(n_series=27)
    cell = make_cell(SPEC, 70.0, eol=True)
    i1 = sd.section_current(cell, 49.98)
    i2 = sd.section_current(cell, 50.02)
    assert 0 < i1 - i2 < 1e-3 * i1


# -------------------------------------------------------------- regulator

def _sim(i_s, n, c, esr, h, s, i_load):
    v0 = 50.0 - s * (n - 1) / 2
    i_sec = np.full(n, i_s)
    th = steady_state_theory(i_sec, i_load, c, esr, h, s, v0)
    sim = S3R(i_sec, c, esr, v0, h, s, vc=th["v_mean"], i_load=i_load)
    sim.shunted[: th["j"] + 1] = True
    return sim, th


@pytest.mark.parametrize("esr", [0.0, 0.004, 0.012])
@pytest.mark.parametrize("x_frac", [0.03, 0.2, 0.5, 0.81, 0.97])
def test_frequency_duty_and_mean_voltage(esr, x_frac):
    i_s, n, c, h, s = 2.1, 12, 3.4e-3, 0.040, 0.050
    sim, th = _sim(i_s, n, c, esr, h, s, (5 + x_frac) * i_s)
    sim.run(12.0 / th["f"])
    period, on_frac, v_mean = measure_cycle(sim, th["j"])
    assert 1.0 / period == pytest.approx(th["f"], rel=1e-9)
    assert on_frac == pytest.approx(x_frac, rel=1e-9)
    assert v_mean == pytest.approx(th["v_mean"], abs=1e-9)
    # Only the active section ever switches, and V_bus never leaves its band.
    switched = {k for _, k, a, _ in sim.log if a != 0}
    assert switched == {th["j"]}
    vb = [e[3] for e in sim.log]
    assert min(vb) >= sim.lo[th["j"]] - 1e-12
    assert max(vb) <= sim.hi[th["j"]] + 1e-12


def test_worst_case_frequency_at_half_section():
    i_s, c, h, esr = 2.125, 3.4e-3, 0.040, 0.004
    xs = np.linspace(0.01, 0.99, 99) * i_s
    f = [steady_state_theory([i_s] * 4, 1 * i_s + x, c, esr, h, 0.05, 49.9)["f"]
         for x in xs]
    assert xs[int(np.argmax(f))] == pytest.approx(i_s / 2)
    assert max(f) == pytest.approx(i_s / (4 * c * (h - esr * i_s)), rel=1e-12)


def test_static_regulation_staircase():
    # Mean bus voltage drops by one window spacing per section of load:
    # sunlight output resistance s / Is.
    i_s, n, c, esr, h, s = 2.0, 10, 2e-3, 0.003, 0.040, 0.050
    loads = (np.arange(2, 8) + 0.5) * i_s
    v = []
    for il in loads:
        sim, th = _sim(i_s, n, c, esr, h, s, il)
        sim.run(10.0 / th["f"])
        v.append(measure_cycle(sim, th["j"])[2])
    slope = np.polyfit(loads, v, 1)[0]
    assert slope == pytest.approx(-s / i_s, rel=1e-9)


def test_load_step_extremes_hit_window_edges():
    i_s, n, c, esr, h, s = 2.1, 30, 3.4e-3, 0.004, 0.040, 0.050
    il0 = 10.3 * i_s
    sim, th0 = _sim(i_s, n, c, esr, h, s, il0)
    th1 = steady_state_theory(sim.i_sec, il0 + 20.0, c, esr, h, s, sim.v0)
    sim.run(8e-3, load_steps=[(1e-3, il0 + 20.0), (4e-3, il0)])
    t = np.array([e[0] for e in sim.log])
    vb = np.array([e[3] for e in sim.log])
    assert vb[(t > 1e-3) & (t < 4e-3)].min() == pytest.approx(sim.lo[th1["j"]], abs=1e-12)
    assert vb[t > 4e-3].max() == pytest.approx(sim.hi[th0["j"]], abs=1e-12)
    # back to the original limit cycle
    period, _, _ = measure_cycle(sim, th0["j"])
    assert 1 / period == pytest.approx(th0["f"], rel=1e-9)


def test_charge_balance_over_transient():
    # C (Vc_end - Vc_start) = integral of capacitor current.
    i_s, n, c, esr, h, s = 1.9, 16, 2.5e-3, 0.005, 0.040, 0.050
    sim, _ = _sim(i_s, n, c, esr, h, s, 6.4 * i_s)
    vc0 = sim.vc
    sim.run(5e-3, load_steps=[(0.7e-3, 12.2 * i_s), (2.9e-3, 3.1 * i_s)])
    tr = sim.trace
    q = sum(ia * (tb - ta) for (ta, _, ia), (tb, _, _) in zip(tr[:-1], tr[1:]))
    assert c * (sim.vc - vc0) == pytest.approx(q, abs=1e-12)


def test_against_fixed_step_integration():
    # Independent check: brute-force Euler with sampled comparators.
    i_s, n, c, esr, h, s = 2.0, 6, 1e-3, 0.004, 0.040, 0.050
    il = 2.37 * i_s
    sim, th = _sim(i_s, n, c, esr, h, s, il)
    lo, hi = sim.lo.copy(), sim.hi.copy()
    shunted = sim.shunted.copy()
    vc, dt = sim.vc, 2e-8
    shunt_times = []
    for k in range(int(4e-3 / dt)):
        ic = i_s * np.count_nonzero(~shunted) - il
        vb = vc + esr * ic
        up = ~shunted & (vb >= hi)
        dn = shunted & (vb <= lo)
        if up.any():
            shunted[np.flatnonzero(up)[0]] = True
            if np.flatnonzero(up)[0] == th["j"]:
                shunt_times.append(k * dt)
            continue
        if dn.any():
            shunted[np.flatnonzero(dn)[-1]] = False
            continue
        vc += ic / c * dt
    f_euler = (len(shunt_times) - 2) / (shunt_times[-1] - shunt_times[1])
    assert f_euler == pytest.approx(th["f"], rel=2e-3)


def test_rejects_esr_chatter_and_overlap():
    with pytest.raises(ValueError):
        S3R(np.full(4, 2.0), 1e-3, 0.025, 49.9, 0.040, 0.050)   # ESR Is > h
    with pytest.raises(ValueError):
        S3R(np.full(4, 2.0), 1e-3, 0.0, 49.9, 0.040, 0.035)     # s < h
