# Astron Built this AP
# Project 052 | Thermal analysis / Atmospheric entry | Ballistic and lifting Earth entry: Allen-Eggers, equilibrium glide and Sutton-Graves heating
# Author: Aditya (aditya-astron)
# Date: 2026-11-21 IST
"""Run with: python -m pytest -q"""
import numpy as np
import pytest
from scipy.integrate import quad

from entry import (G0, K_SG, R_E, Planet, Vehicle, ae_heat_load, ae_peak_decel,
                   ae_peak_heating, ae_speed_at_zero_density, ae_velocity,
                   eq_glide_altitude, eq_glide_gamma, eq_glide_range, fly)

AE = Planet(gravity=False, curvature=False)
FULL = Planet()
H0 = 120e3


@pytest.mark.parametrize("gdeg", [-3.0, -12.0, -60.0])
@pytest.mark.parametrize("Ve", [7.8e3, 11.0e3])
def test_allen_eggers_matches_integrator(gdeg, Ve):
    veh = Vehicle(300.0, 0.0, 0.5)
    g = np.radians(gdeg)
    rho_e = AE.rho(H0)
    tr = fly(veh, AE, Ve, g, h0=H0, V_end=1e-3)
    V0 = ae_speed_at_zero_density(Ve, veh.beta, g, rho_e)

    # velocity-density relation along the whole path
    assert np.allclose(tr.V, ae_velocity(AE.rho(tr.h), Ve, veh.beta, g, rho_e), rtol=1e-9)
    # gamma frozen without gravity and curvature
    assert np.allclose(tr.gamma, g, atol=1e-12)

    a, ha = ae_peak_decel(V0, veh.beta, g)
    an, han, Van = tr.peak("decel")
    assert an == pytest.approx(a, rel=1e-8)
    assert han == pytest.approx(ha, abs=1.0)
    assert Van == pytest.approx(V0 / np.sqrt(np.e), rel=1e-6)

    q, hq, Vq = ae_peak_heating(V0, veh.beta, g, veh.nose_radius)
    qn, hqn, Vqn = tr.peak("qdot")
    assert qn == pytest.approx(q, rel=1e-8)
    assert hqn == pytest.approx(hq, abs=1.0)
    assert Vqn == pytest.approx(Vq, rel=1e-6)

    Q = ae_heat_load(V0, veh.beta, g, veh.nose_radius, rho_e=rho_e)
    assert tr.heat_load() == pytest.approx(Q, rel=1e-8)


def test_heat_load_closed_form_by_quadrature():
    # independent check of the erf form: integrate k sqrt(rho/rn) V^2 H / (rho s) d rho
    Ve, beta, g, rn = 7.8e3, 300.0, np.radians(-8.0), 0.5
    s, H = np.sin(-g), AE.scale_height
    f = lambda rho: K_SG * np.sqrt(rho / rn) * ae_velocity(rho, Ve, beta, g) ** 2 * H / (rho * s)
    num = quad(f, 0.0, 1e-4, epsabs=0, epsrel=1e-12)[0] + quad(f, 1e-4, 1.225, epsabs=0, epsrel=1e-12)[0]
    assert ae_heat_load(Ve, beta, g, rn) == pytest.approx(num, rel=1e-9)


def test_hand_numbers():
    # Ve = 7.8 km/s, 20 deg: a = V^2 sin g / (2 e H) = 6.084e7 * 0.34202 / 39144.5 = 531.6 m/s^2
    a, _ = ae_peak_decel(7.8e3, 300.0, np.radians(-20.0))
    assert a == pytest.approx(531.6, rel=1e-3)
    assert a / G0 == pytest.approx(54.21, abs=0.01)
    # beta = 300, s = 0.34202 -> rho* = 0.014251 kg/m^3 -> h = 7.2 ln(1.225/0.014251) = 32.06 km
    _, h = ae_peak_decel(7.8e3, 300.0, np.radians(-20.0))
    assert h == pytest.approx(32.06e3, abs=20)


def test_steep_entry_gravity_bias():
    # Vertical full-model entry: gravity adds about 2 g dh / Ve^2 to V^2 by peak g,
    # so the closed form should be 2.5-4 % low, not more.
    veh = Vehicle(300.0, 0.0, 0.5)
    tr = fly(veh, FULL, 7.8e3, np.radians(-90.0), h0=H0)
    an, h_pk, _ = tr.peak("decel")
    a, _ = ae_peak_decel(7.8e3, veh.beta, np.radians(-90.0))
    g_avg = FULL.g(0.5 * (H0 + h_pk))
    expected = 2 * g_avg * (H0 - h_pk) / 7.8e3 ** 2
    assert an / a - 1 == pytest.approx(expected, abs=0.006)


def test_peak_g_independent_of_beta_at_circular_speed():
    peaks = [fly(Vehicle(b, 0.0, 0.5), FULL, 7.8e3, np.radians(-6.0), h0=H0).peak("decel")[0]
             for b in (100.0, 300.0, 1000.0)]
    assert max(peaks) / min(peaks) - 1 < 0.01


def test_shallow_lunar_return_skips():
    tr = fly(Vehicle(300.0, 0.0, 0.5), FULL, 11.0e3, np.radians(-2.0), h0=H0)
    assert tr.skipped
    tr = fly(Vehicle(300.0, 0.0, 0.5), FULL, 11.0e3, np.radians(-6.0), h0=H0)
    assert not tr.skipped


def test_equilibrium_glide_state_and_range():
    veh = Vehicle(400.0, 2.0, 0.3)
    V0 = 7.5e3
    h0 = eq_glide_altitude(V0, veh, FULL)
    r = R_E + h0
    # lift balances weight less centrifugal relief at the computed altitude
    lift = FULL.rho(h0) * V0 ** 2 / (2 * veh.beta) * veh.lift_to_drag
    assert lift == pytest.approx(FULL.g(h0) - V0 ** 2 / r, rel=1e-10)
    g0 = eq_glide_gamma(V0, veh.lift_to_drag, h0, FULL)
    tr = fly(veh, FULL, V0, g0, h0=h0, V_end=1.5e3)
    s_cf = eq_glide_range(V0, 1.5e3, veh.lift_to_drag, r, FULL.g(h0)) * R_E / r
    assert tr.s[-1] / s_cf - 1 == pytest.approx(0.0, abs=0.03)
    # started on the equilibrium path, so altitude stays within a few km of it early on
    i = np.searchsorted(-tr.V, -6.0e3)
    h_eq = eq_glide_altitude(tr.V[i], veh, FULL)
    assert abs(tr.h[i] - h_eq) < 3e3
