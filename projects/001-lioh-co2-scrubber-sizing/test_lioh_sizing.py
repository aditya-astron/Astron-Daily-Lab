# Astron Built this AP
# Project 001 | ECLSS | LiOH CO2 Scrubber Sizing for a 4-Crew 30-Day Mission
# Author: Aditya (aditya-astron)
# Date: 2026-10-30 IST
"""Checks against hand calculations and published reference values."""

import math

import numpy as np
import pytest

from lioh_sizing import (Canister, Mission, kg_m3_to_mmhg, lioh_per_kg_co2,
                         pa_to_kg_m3, reaction_heat_per_kg_co2, simulate,
                         static_budget, steady_state_ppco2, water_per_kg_co2)


def test_stoichiometry():
    # 2 * 23.946 / 44.009 = 1.0882; commonly quoted as ~1.09 kg LiOH / kg CO2
    assert lioh_per_kg_co2() == pytest.approx(1.0882, abs=1e-4)
    assert water_per_kg_co2() == pytest.approx(0.4093, abs=1e-4)


def test_reaction_enthalpy():
    # dH = (-1216.0 - 241.8) - (2*-484.9 - 393.5) = -94.5 kJ/mol CO2
    assert reaction_heat_per_kg_co2() == pytest.approx(94.5e3 / 44.009e-3, rel=1e-6)
    assert reaction_heat_per_kg_co2(water_liquid=True) > reaction_heat_per_kg_co2()


def test_unit_round_trip():
    c = pa_to_kg_m3(400.0, 295.15)
    assert kg_m3_to_mmhg(c, 295.15) == pytest.approx(400.0 / 133.322368)


def test_static_budget():
    b = static_budget(Mission())
    assert b["co2_kg"] == pytest.approx(4 * 1.04 * 30)
    assert b["lioh_kg"] == pytest.approx(124.8 * 1.0882, rel=1e-3)


def test_steady_state_hand_calc():
    # C = G/(eta Q) = 4.815e-5 / (0.98 * 0.01) kg/m^3 -> p = C R T / M
    m, c = Mission(), Canister()
    conc = (4.16 / 86400) / (0.98 * 0.01)
    p_mmhg = conc * 8.314462618 * 295.15 / 44.009e-3 / 133.322368
    assert steady_state_ppco2(m, c) == pytest.approx(p_mmhg, rel=1e-9)
    assert 2.0 < p_mmhg < 2.1


def test_efficiency_curve():
    c = Canister()
    assert c.efficiency(0.0) == pytest.approx(c.eta_fresh)
    assert c.efficiency(c.u50) == pytest.approx(0.5 * c.eta_fresh, rel=1e-6)
    assert c.efficiency(1.0) < 0.02


def test_fresh_bed_reaches_steady_state():
    # Huge canister: never swaps, cabin settles at G/(eta Q) with
    # time constant V/(eta Q) ~ 34 min.
    m = Mission(days=1.0)
    c = Canister(lioh_mass=1000.0)
    sim = simulate(m, c)
    assert sim.canisters_used == 1
    assert sim.ppco2[-1] == pytest.approx(steady_state_ppco2(m, c), rel=1e-4)


def test_mass_conservation_and_limit():
    m, c = Mission(), Canister()
    sim = simulate(m, c)
    err = sim.co2_generated - sim.co2_absorbed - sim.co2_cabin_delta
    assert abs(err) < 1e-6
    assert sim.ppco2.max() <= m.ppco2_limit_mmhg + 1e-6


def test_canister_count_matches_capacity_estimate():
    m, c = Mission(), Canister()
    sim = simulate(m, c)
    u = np.mean(sim.swap_utilisation)
    # Each swapped canister held u * capacity; total must cover the mission.
    expected = math.ceil(sim.co2_absorbed / (u * c.co2_capacity) - 1e-9)
    assert abs(sim.canisters_used - expected) <= 1
    # Canister life from hand calc: u * cap / G
    life_h = u * c.co2_capacity / (m.co2_rate * 3600)
    assert np.mean(np.diff([0.0] + sim.swap_times)) / 3600 == pytest.approx(life_h, rel=0.02)
