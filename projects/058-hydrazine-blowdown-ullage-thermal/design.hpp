// Astron Built this AP
// Project 058 | Rocket propulsion | Hydrazine blowdown system sizing with ullage-gas thermodynamics
// Author: Aditya (aditya-astron)
// Date: 2026-11-23 IST
//
// Mission requirements and the isothermal blowdown sizing that conventional
// first-pass design uses, plus the bisection that re-sizes the loaded pressure
// once the ullage thermodynamics are flown.

#pragma once

#include "blowdown.hpp"

namespace bd {

struct Requirements {
    double m_sc0 = 180.0;          // wet spacecraft mass, kg
    double dv = 120.0;             // m/s
    double residual = 0.03;        // unexpellable fraction of the load
    double p_min = 5.5 * kBar;     // lowest feed pressure the thrusters are rated for
    double p_max = 24.0 * kBar;    // highest feed pressure the thrusters accept
    double burn = 900.0;           // s per manoeuvre
    double coast = 5400.0;         // s between manoeuvres (about one LEO orbit)
    double thrust_nom = 1.0;       // N per thruster at 22 bar feed
};

inline Thruster baseline_thruster(double thrust_nom = 1.0) {
    // Nominal point: 22 bar feed with 11 bar in the chamber.
    return Thruster::from_nominal(thrust_nom, 22.0 * kBar, 11.0 * kBar);
}

struct Budget {
    double m_usable, m_load, v_usable, v_load;
};

inline Budget propellant_budget(const Requirements& r, const Thruster& thr, double rho_p) {
    Budget b{};
    b.m_usable = r.m_sc0 * (1.0 - std::exp(-r.dv / (thr.isp() * kG0)));
    b.m_load = b.m_usable / (1.0 - r.residual);
    b.v_usable = b.m_usable / rho_p;
    b.v_load = b.m_load / rho_p;
    return b;
}

// Isothermal sizing: Pend = p_min exactly at depletion for blowdown ratio B.
inline System isothermal_design(const Requirements& r, double blowdown, double h) {
    System s;
    s.thr = baseline_thruster(r.thrust_nom);
    const Budget b = propellant_budget(r, s.thr, s.rho_p);
    s.v_usable = b.v_usable;
    s.vg0 = b.v_usable / (blowdown - 1.0);
    s.vt = b.v_load + s.vg0;
    s.p0 = blowdown * r.p_min;
    s.m_sc0 = r.m_sc0;
    s.h = h;
    return s;
}

// Same tank and load, more pressurant: find P0 such that the lowest feed
// pressure seen anywhere in the manoeuvre plan equals p_min. Returns a
// negative value when even p_max cannot hold the floor.
inline double resize_p0(System s, const Requirements& r) {
    auto floor_at = [&](double p0) {
        s.p0 = p0;
        return fly_mission(s, r.burn, r.coast).p_min;
    };
    double lo = r.p_min, hi = r.p_max;
    if (floor_at(hi) < r.p_min) return -1.0;
    for (int i = 0; i < 50; ++i) {
        const double mid = 0.5 * (lo + hi);
        if (floor_at(mid) < r.p_min) lo = mid; else hi = mid;
    }
    return hi;
}

// Alternative fix when P0 is capped: keep the loaded pressure, grow the ullage
// (and the tank) until the floor holds.
inline double resize_vg0(System s, const Requirements& r) {
    const double vl = s.vt - s.vg0;
    auto floor_at = [&](double vg0) {
        s.vg0 = vg0;
        s.vt = vl + vg0;
        return fly_mission(s, r.burn, r.coast).p_min;
    };
    double lo = s.vg0, hi = 2.0 * s.vg0;
    while (floor_at(hi) < r.p_min) { lo = hi; hi *= 2.0; }
    for (int i = 0; i < 50; ++i) {
        const double mid = 0.5 * (lo + hi);
        if (floor_at(mid) < r.p_min) lo = mid; else hi = mid;
    }
    return hi;
}

}  // namespace bd
