# Astron Built this AP
# Project 052 | Thermal analysis / Atmospheric entry | Ballistic and lifting Earth entry: Allen-Eggers, equilibrium glide and Sutton-Graves heating
# Author: Aditya (aditya-astron)
# Date: 2026-11-21 IST
"""
Entry study, all starting at the 120 km interface:

  1. Allen-Eggers closed forms vs. the integrator with gravity and curvature
     switched off (should agree to solver tolerance).
  2. The same closed forms vs. the full model: where the straight-line,
     no-gravity assumption stops being good enough.
  3. Ballistic trade over entry angle and beta for LEO return (7.8 km/s) and
     lunar return (11.0 km/s): peak g, peak heat flux, heat load, pulse width.
  4. Equilibrium glide: closed-form range vs. the full model from the
     equilibrium state, for L/D = 0.3, 1 and 2.
"""
import numpy as np

from entry import (G0, R_E, Planet, Vehicle, ae_heat_load, ae_peak_decel,
                   ae_peak_heating, ae_speed_at_zero_density, eq_glide_altitude,
                   eq_glide_gamma, eq_glide_range, fly)

H0 = 120e3
FULL = Planet()
AE_WORLD = Planet(gravity=False, curvature=False)


def pulse_width(tr, frac=0.5):
    """Time the heat flux spends above frac x peak, s."""
    above = tr.t[tr.qdot >= frac * tr.qdot.max()]
    return above[-1] - above[0]


def section_closed_form_check():
    print("1. Allen-Eggers closed forms vs. integrator (no gravity, flat, beta = 300 kg/m^2, rn = 0.5 m)")
    print(f"   {'gamma':>6} {'Ve km/s':>8} {'peak g':>8} {'rel err':>9} {'h_gmax km':>10} "
          f"{'q MW/m2':>8} {'rel err':>9} {'Q MJ/m2':>8} {'rel err':>9}")
    veh = Vehicle(300.0, 0.0, 0.5)
    rho_e = AE_WORLD.rho(H0)
    for Ve in (7.8e3, 11.0e3):
        for gdeg in (-3.0, -10.0, -45.0):
            g = np.radians(gdeg)
            tr = fly(veh, AE_WORLD, Ve, g, h0=H0, V_end=1e-3)
            V0 = ae_speed_at_zero_density(Ve, veh.beta, g, rho_e)
            a, ha = ae_peak_decel(V0, veh.beta, g)
            q, _, _ = ae_peak_heating(V0, veh.beta, g, veh.nose_radius)
            Q = ae_heat_load(V0, veh.beta, g, veh.nose_radius, rho_e=rho_e)
            an, _, _ = tr.peak("decel")
            qn, _, _ = tr.peak("qdot")
            Qn = tr.heat_load()
            print(f"   {gdeg:6.1f} {Ve/1e3:8.1f} {an/G0:8.2f} {an/a-1:9.1e} {ha/1e3:10.2f} "
                  f"{qn/1e6:8.3f} {qn/q-1:9.1e} {Qn/1e6:8.1f} {Qn/Q-1:9.1e}")
    print()


def section_full_vs_closed_form():
    print("2. Full model (gravity + curvature) vs. Allen-Eggers, Ve = 7.8 km/s, beta = 300, rn = 0.5")
    print(f"   {'gamma':>6} {'g full':>7} {'g AE':>7} {'err':>7} {'q full':>7} {'q AE':>7} {'err':>7} "
          f"{'h_q full':>8} {'h_q AE':>7}")
    veh = Vehicle(300.0, 0.0, 0.5)
    for gdeg in (-1.5, -2.5, -5.0, -10.0, -20.0, -45.0, -90.0):
        g = np.radians(gdeg)
        tr = fly(veh, FULL, 7.8e3, g, h0=H0)
        a, _ = ae_peak_decel(7.8e3, veh.beta, g)
        q, hq, _ = ae_peak_heating(7.8e3, veh.beta, g, veh.nose_radius)
        an, _, _ = tr.peak("decel")
        qn, hqn, _ = tr.peak("qdot")
        print(f"   {gdeg:6.1f} {an/G0:7.2f} {a/G0:7.2f} {an/a-1:+7.1%} {qn/1e6:7.3f} {q/1e6:7.3f} "
              f"{qn/q-1:+7.1%} {hqn/1e3:8.1f} {hq/1e3:7.1f}")
    print("   Shallow entries: gravity and the falling-away surface steepen the path, so the vehicle")
    print("   meets the dense layer steeper than gamma_e and AE badly underpredicts the load. Steep")
    print("   entries: the closed forms run about 3 % low in g and 4 % low in flux, the kinetic")
    print("   energy gravity adds on the way down (2 g dh / Ve^2 is about 3 % over 90 km).")
    print()


def section_trade():
    print("3. Ballistic trade (full model), rn = 0.5 m")
    print(f"   {'Ve':>5} {'beta':>5} {'gamma':>6} {'peak g':>7} {'q MW/m2':>8} {'Q MJ/m2':>8} "
          f"{'q>50% s':>8} {'h_qmax km':>9}")
    for Ve in (7.8e3, 11.0e3):
        for beta in (100.0, 300.0, 1000.0):
            for gdeg in (-2.0, -6.0, -20.0):
                veh = Vehicle(beta, 0.0, 0.5)
                tr = fly(veh, FULL, Ve, np.radians(gdeg), h0=H0)
                if tr.skipped:
                    print(f"   {Ve/1e3:5.1f} {beta:5.0f} {gdeg:6.1f}   skips out: min h "
                          f"{tr.h.min()/1e3:.1f} km, exit speed {tr.V[-1]/1e3:.2f} km/s")
                    continue
                an, _, _ = tr.peak("decel")
                qn, hq, _ = tr.peak("qdot")
                print(f"   {Ve/1e3:5.1f} {beta:5.0f} {gdeg:6.1f} {an/G0:7.2f} {qn/1e6:8.3f} "
                      f"{tr.heat_load()/1e6:8.1f} {pulse_width(tr):8.1f} {hq/1e3:9.1f}")
    print("   At 7.8 km/s peak g depends on Ve and gamma only (beta just moves where it happens),")
    print("   peak flux grows as sqrt(beta sin gamma) and heat load as sqrt(beta / sin gamma): steep")
    print("   entries give short intense pulses, shallow ones long soaks into the TPS. At 11 km/s the")
    print("   speed is supercircular: a 2 deg ballistic entry skips straight out, and at 6 deg a")
    print("   high-beta body penetrates deeper before it slows, trading peak g for flux and load.")
    print()


def section_glide():
    print("4. Equilibrium glide from 7.5 km/s to 1.5 km/s, beta = 400 kg/m^2")
    print(f"   {'L/D':>4} {'h0 km':>6} {'gamma0 deg':>10} {'range km':>9} {'closed km':>9} {'err':>7} "
          f"{'peak g':>7} {'eq g':>6}")
    V0, Vf = 7.5e3, 1.5e3
    for LD in (0.3, 1.0, 2.0):
        veh = Vehicle(400.0, LD, 0.3)
        h0 = eq_glide_altitude(V0, veh, FULL)
        g0 = eq_glide_gamma(V0, LD, h0, FULL)
        tr = fly(veh, FULL, V0, g0, h0=h0, V_end=Vf)
        r = R_E + h0
        s_cf = eq_glide_range(V0, Vf, LD, r, FULL.g(h0)) * R_E / r
        an, _, _ = tr.peak("decel")
        a_eq_end = FULL.g(h0) * (1 - Vf ** 2 / (FULL.g(h0) * r)) / LD * np.hypot(1, LD)
        print(f"   {LD:4.1f} {h0/1e3:6.1f} {np.degrees(g0):10.3f} {tr.s[-1]/1e3:9.0f} {s_cf/1e3:9.0f} "
              f"{tr.s[-1]/s_cf-1:+7.1%} {an/G0:7.2f} {a_eq_end/G0:6.2f}")
    print("   'eq g' is the total aero load at Vf on the equilibrium path, its largest value; the")
    print("   flown peaks stay under it. The closed form drops the g sin(gamma) energy gained while")
    print("   sinking and the phugoid the start state excites; the range error is 1-2 % at L/D of")
    print("   1-2 and 8 % at a capsule's 0.3, so treat it as a guide only for real gliders.")
    print()


if __name__ == "__main__":
    section_closed_form_check()
    section_full_vs_closed_form()
    section_trade()
    section_glide()
