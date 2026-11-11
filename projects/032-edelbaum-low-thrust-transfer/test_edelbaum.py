# Astron Built this AP
# Project 032 | Orbital mechanics | Edelbaum low-thrust LEO-to-GEO transfer with plane change
# Author: Aditya (aditya-astron)
# Date: 2026-11-11 IST
"""Verification of the Edelbaum transfer model. Run: python3 -m pytest -q"""

import math

import pytest

from edelbaum import (MU, R_EARTH, R_GEO, averaged_check, circ_speed,
                      direct_optimum, edelbaum, hohmann_plane_change,
                      propagate_cartesian, sequential_lowthrust, v_of_dv)

R_LEO = R_EARTH + 200.0
V_LEO, V_GEO = circ_speed(R_LEO), circ_speed(R_GEO)
DI = math.radians(28.5)


def test_coplanar_limit_is_speed_difference():
    dv, b0, b1 = edelbaum(V_LEO, V_GEO, 0.0)
    assert dv == pytest.approx(V_LEO - V_GEO, rel=1e-14)
    assert b0 == pytest.approx(0.0, abs=1e-15)


def test_equal_radius_beats_pure_plane_change():
    # Same orbit, 10 deg change: Edelbaum climbs where turning is cheaper,
    # so it beats thrusting purely out of plane (pi/2 v di).
    di = math.radians(10.0)
    dv, b0, _ = edelbaum(V_LEO, V_LEO, di)
    assert dv == pytest.approx(2.0 * V_LEO * math.sin(math.pi * di / 4.0), rel=1e-13)
    assert dv < 0.5 * math.pi * V_LEO * di
    assert b0 < 0.5 * math.pi


def test_leo_geo_hand_value():
    # sqrt(v0^2 - 2 v0 v1 cos(pi*di/2) + v1^2) evaluated by hand:
    # v0 = 7.78426, v1 = 3.07466, cos(0.78147) = 0.70991  ->  6.0053 km/s
    assert V_LEO == pytest.approx(7.78426, abs=1e-5)
    assert V_GEO == pytest.approx(3.07466, abs=1e-5)
    dv, b0, b1 = edelbaum(V_LEO, V_GEO, DI)
    assert dv == pytest.approx(6.0053, abs=1e-4)
    assert math.degrees(b0) == pytest.approx(21.135, abs=1e-3)
    # End state from the closed-form schedule matches the target orbit.
    assert v_of_dv(V_LEO, b0, dv) == pytest.approx(V_GEO, rel=1e-12)
    assert 2 / math.pi * (b1 - b0) == pytest.approx(DI, rel=1e-12)


def test_averaged_ode_reproduces_closed_form():
    dv, b0, _ = edelbaum(V_LEO, V_GEO, DI)
    v_end, di_end = averaged_check(V_LEO, b0, dv)
    assert v_end == pytest.approx(V_GEO, abs=1e-9)
    assert di_end == pytest.approx(DI, abs=1e-10)


@pytest.mark.parametrize("v0,v1,di_deg", [
    (V_LEO, V_GEO, 28.5),
    (V_LEO, circ_speed(26560.0), 55.0),   # LEO to GPS-like MEO
    (circ_speed(9000.0), circ_speed(7000.0), 20.0),  # lowering, beta > 90 deg
])
def test_direct_optimisation_cannot_beat_edelbaum(v0, v1, di_deg):
    di = math.radians(di_deg)
    dv_e, _, _ = edelbaum(v0, v1, di)
    dv_d, _, ok = direct_optimum(v0, v1, di, n_seg=30)
    assert ok
    assert dv_d >= dv_e * (1 - 1e-9)
    assert dv_d <= dv_e * (1 + 5e-4)


def test_sequential_is_worse():
    dv, _, _ = edelbaum(V_LEO, V_GEO, DI)
    assert sequential_lowthrust(V_LEO, V_GEO, DI) > dv


def test_hohmann_coplanar_and_split():
    a = 0.5 * (R_LEO + R_GEO)
    ref = (math.sqrt(MU / R_LEO) * (math.sqrt(R_GEO / a) - 1.0)
           + math.sqrt(MU / R_GEO) * (1.0 - math.sqrt(R_LEO / a)))
    dv0, alpha0 = hohmann_plane_change(R_LEO, R_GEO, 0.0)
    assert dv0 == pytest.approx(ref, rel=1e-10)
    # Hand: 7.78426 * 0.31533 + 3.07466 * 0.48047 = 2.4546 + 1.4773
    assert ref == pytest.approx(3.9319, abs=2e-4)
    dv, alpha = hohmann_plane_change(R_LEO, R_GEO, DI)
    # Most of the plane change belongs at apogee, a couple of degrees at perigee.
    assert 1.5 < math.degrees(alpha) < 3.0
    # Doing the whole plane change at apogee costs more than the split.
    va = math.sqrt(MU * (2.0 / R_GEO - 1.0 / a))
    vp = math.sqrt(MU * (2.0 / R_LEO - 1.0 / a))
    all_apogee = (vp - V_LEO) + math.sqrt(va**2 + V_GEO**2 - 2 * va * V_GEO * math.cos(DI))
    assert dv < all_apogee
    edel, _, _ = edelbaum(V_LEO, V_GEO, DI)
    assert dv < edel


@pytest.mark.parametrize("r0,r1,i0_deg,i1_deg", [
    (7000.0, 10000.0, 10.0, 0.0),   # raise and lower inclination
    (10000.0, 7000.0, 5.0, 15.0),   # lower orbit, raise inclination
])
def test_cartesian_flight_reaches_target(r0, r1, i0_deg, i1_deg):
    res = propagate_cartesian(r0, r1, math.radians(i0_deg), math.radians(i1_deg),
                              thrust_n=1.0, m0_kg=1000.0, isp_s=1600.0)
    assert res["status"] == 1                     # stopped on the dV event
    assert res["dv_flown"] == pytest.approx(res["dv_target"], rel=1e-9)
    assert res["t_end"] == pytest.approx(res["t_rocket"], rel=1e-6)
    assert res["a"] == pytest.approx(r1, rel=1e-3)
    assert math.degrees(res["i"]) == pytest.approx(i1_deg, abs=0.02)
    assert res["e"] < 5e-3
