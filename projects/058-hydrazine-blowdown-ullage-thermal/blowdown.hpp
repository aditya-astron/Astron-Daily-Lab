// Astron Built this AP
// Project 058 | Rocket propulsion | Hydrazine blowdown system sizing with ullage-gas thermodynamics
// Author: Aditya (aditya-astron)
// Date: 2026-11-23 IST
//
// Models for a monopropellant blowdown feed: a lumped thruster flow model, a
// spherical diaphragm tank with an N2 ullage, and the ullage energy balance
// that sits between the isothermal and adiabatic limits.

#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

namespace bd {

constexpr double kPi = 3.14159265358979323846;
constexpr double kG0 = 9.80665;
constexpr double kBar = 1.0e5;

// ---------------------------------------------------------------------------
// Thruster. Catalyst-bed monoprop, all feed/valve/injector/bed losses lumped
// into one orifice:   mdot = K sqrt(Pf - Pc)
// and the nozzle throat chokes:  mdot = Pc At / c*.
// Equating them gives a quadratic in Pc:
//     a Pc^2 + K^2 Pc - K^2 Pf = 0,   a = (At / c*)^2
// Thrust F = Cf Pc At, so Isp = c* Cf / g0 is pressure independent here. Real
// thrusters lose a few seconds of Isp at the bottom of blowdown; that changes
// propellant load by well under 1 % and is left out.
// ---------------------------------------------------------------------------
struct Thruster {
    double cstar = 1320.0;  // m/s, hydrazine with ~40 % NH3 dissociation
    double cf = 1.66;       // vacuum thrust coefficient, area ratio ~ 50-100
    double at = 0.0;        // m^2 throat area
    double k = 0.0;         // kg/s/sqrt(Pa) lumped feed conductance

    // Size the throat and feed conductance from one nominal operating point.
    static Thruster from_nominal(double thrust_n, double pf_pa, double pc_pa) {
        Thruster t;
        t.at = thrust_n / (t.cf * pc_pa);
        const double mdot = pc_pa * t.at / t.cstar;
        t.k = mdot / std::sqrt(pf_pa - pc_pa);
        return t;
    }

    double isp() const { return cstar * cf / kG0; }

    // Stable root of the quadratic (written to avoid cancellation at small Pf).
    double chamber_pressure(double pf) const {
        if (pf <= 0.0) return 0.0;
        const double a = (at / cstar) * (at / cstar);
        const double k2 = k * k;
        return 2.0 * k2 * pf / (k2 + std::sqrt(k2 * k2 + 4.0 * a * k2 * pf));
    }
    double mdot(double pf) const { return chamber_pressure(pf) * at / cstar; }
    double thrust(double pf) const { return cf * chamber_pressure(pf) * at; }
};

// ---------------------------------------------------------------------------
// Pressurant: nitrogen, calorically perfect over 230-300 K.
// ---------------------------------------------------------------------------
struct Gas {
    double r = 296.8;  // J/(kg K)
    double gamma = 1.4;
    double cv() const { return r / (gamma - 1.0); }
};

// ---------------------------------------------------------------------------
// Blowdown sizing. Let Vl = loaded propellant volume, Vu = usable (expelled)
// volume, Vg0 = initial ullage, B = P0/Pend the blowdown ratio, and assume a
// polytropic expansion P V^n = const over the whole mission. Then
//     Vg0 = Vu / (B^(1/n) - 1),   Vt = Vl + Vg0
// For a pressure-sized thin sphere tank mass ~ P0 Vt, i.e. proportional to
//     g(B) = B (Vl + Vu / (B^(1/n) - 1))      (Pend fixed by the thruster)
// Closed forms of the optimum:
//     n = 1:           B* = 1 + sqrt(Vu/Vl)    (= 2 with no residual)
//     Vl = Vu, any n:  B* = ((n+1)/n)^n        (-> e as n -> inf)
// ---------------------------------------------------------------------------
inline double pv_metric(double b, double vl, double vu, double n) {
    return b * (vl + vu / (std::pow(b, 1.0 / n) - 1.0));
}

inline double golden_min(const std::function<double(double)>& f, double lo, double hi,
                         double tol = 1e-12) {
    const double r = 0.5 * (std::sqrt(5.0) - 1.0);
    double a = lo, b = hi;
    double c = b - r * (b - a), d = a + r * (b - a);
    double fc = f(c), fd = f(d);
    while (b - a > tol * (1.0 + std::fabs(a) + std::fabs(b))) {
        if (fc < fd) {
            b = d; d = c; fd = fc;
            c = b - r * (b - a); fc = f(c);
        } else {
            a = c; c = d; fc = fd;
            d = a + r * (b - a); fd = f(d);
        }
    }
    return 0.5 * (a + b);
}

// Spherical tank shell: membrane stress sizing at burst = sf * MEOP, floored by
// a minimum manufacturable gauge.
struct TankMaterial {
    double rho = 4430.0;      // Ti-6Al-4V, kg/m^3
    double sigma_u = 895e6;   // ultimate tensile, Pa
    double burst_factor = 1.5;
    double min_gauge = 0.6e-3;  // m, machined/spun Ti hemisphere floor
};

struct ShellMass {
    double radius, t_pressure, t_used, mass;
    bool gauge_limited;
};

inline ShellMass shell_mass(double volume, double meop, const TankMaterial& m) {
    ShellMass s{};
    s.radius = std::cbrt(3.0 * volume / (4.0 * kPi));
    s.t_pressure = m.burst_factor * meop * s.radius / (2.0 * m.sigma_u);
    s.t_used = std::max(s.t_pressure, m.min_gauge);
    s.gauge_limited = s.t_used > s.t_pressure;
    s.mass = m.rho * 4.0 * kPi * s.radius * s.radius * s.t_used;
    return s;
}

// ---------------------------------------------------------------------------
// Ullage thermodynamics. Fixed gas mass mg in volume Vg (grows as propellant is
// expelled through the diaphragm). First law for the closed gas system:
//     mg cv dT/dt = -P dVg/dt + h A_g (Tw - T),   P = mg R T / Vg
// h -> 0 is the adiabatic limit (P Vg^gamma = const), h -> inf the isothermal
// one (P Vg = const). Tw is the tank wall/propellant temperature, held fixed:
// ~10 kg of hydrazine plus the shell carry ~40 kJ/K against ~60 J/K of gas.
// The gas-wetted wall area is taken proportional to the ullage fraction,
// A_g = 4 pi r^2 Vg/Vt, which is reasonable for a hemispherical diaphragm.
// ---------------------------------------------------------------------------
struct System {
    Thruster thr;
    Gas gas;
    int n_thrusters = 4;
    double rho_p = 1004.0;  // hydrazine at 20 C, kg/m^3
    double vt = 0.0;        // tank volume
    double vg0 = 0.0;       // initial ullage
    double p0 = 0.0;        // initial pressure at Tw
    double tw = 293.15;
    double h = 5.0;          // W/(m^2 K); <0 selects the isothermal limit
    double v_usable = 0.0;   // propellant volume available for expulsion
    double m_sc0 = 0.0;      // wet spacecraft mass at start

    bool isothermal() const { return h < 0.0; }
    double gas_mass() const { return p0 * vg0 / (gas.r * tw); }
    double radius() const { return std::cbrt(3.0 * vt / (4.0 * kPi)); }
    double wetted_area(double vg) const { return 4.0 * kPi * radius() * radius() * vg / vt; }
};

struct State {
    double t = 0.0;
    double vg = 0.0;
    double temp = 0.0;
    double m_sc = 0.0;
    double dv = 0.0;
};

inline double pressure(const System& s, const State& x) {
    return s.gas_mass() * s.gas.r * x.temp / x.vg;
}

struct Deriv {
    double vg, temp, m_sc, dv;
};

inline Deriv deriv(const System& s, const State& x, bool firing) {
    const double p = pressure(s, x);
    const double mdot = firing ? s.n_thrusters * s.thr.mdot(p) : 0.0;
    const double f = firing ? s.n_thrusters * s.thr.thrust(p) : 0.0;
    Deriv d{};
    d.vg = mdot / s.rho_p;
    if (!s.isothermal()) {
        const double q = s.h * s.wetted_area(x.vg) * (s.tw - x.temp);
        d.temp = (-p * d.vg + q) / (s.gas_mass() * s.gas.cv());
    }
    d.m_sc = -mdot;
    d.dv = f / x.m_sc;
    return d;
}

inline State rk4_step(const System& s, const State& x, double dt, bool firing) {
    auto add = [](const State& a, const Deriv& d, double c) {
        State r = a;
        r.vg += c * d.vg;
        r.temp += c * d.temp;
        r.m_sc += c * d.m_sc;
        r.dv += c * d.dv;
        r.t += c;
        return r;
    };
    const Deriv k1 = deriv(s, x, firing);
    const Deriv k2 = deriv(s, add(x, k1, 0.5 * dt), firing);
    const Deriv k3 = deriv(s, add(x, k2, 0.5 * dt), firing);
    const Deriv k4 = deriv(s, add(x, k3, dt), firing);
    State r = x;
    r.vg += dt / 6.0 * (k1.vg + 2 * k2.vg + 2 * k3.vg + k4.vg);
    r.temp += dt / 6.0 * (k1.temp + 2 * k2.temp + 2 * k3.temp + k4.temp);
    r.m_sc += dt / 6.0 * (k1.m_sc + 2 * k2.m_sc + 2 * k3.m_sc + k4.m_sc);
    r.dv += dt / 6.0 * (k1.dv + 2 * k2.dv + 2 * k3.dv + k4.dv);
    r.t += dt;
    return r;
}

inline State initial_state(const System& s) {
    return State{0.0, s.vg0, s.tw, s.m_sc0, 0.0};
}

struct BurnRecord {
    double t_start, duration, p_start, p_end, t_gas_end, dv;
};

struct MissionResult {
    std::vector<BurnRecord> burns;
    double p_min = 1e30;
    double t_gas_min = 1e30;
    double fire_time = 0.0;
    State final_state;
};

// Fire all thrusters for `burn` seconds, coast `coast` seconds, repeat until the
// usable propellant is gone. The depletion instant is located inside the step
// by bisection on the step length (Vg is monotone in it).
inline MissionResult fly_mission(const System& s, double burn, double coast, double dt = 0.25) {
    MissionResult res;
    State x = initial_state(s);
    const double vg_end = s.vg0 + s.v_usable;
    bool depleted = false;
    while (!depleted) {
        BurnRecord rec{};
        rec.t_start = x.t;
        rec.p_start = pressure(s, x);
        const double dv0 = x.dv;
        double elapsed = 0.0;
        while (elapsed < burn - 1e-12) {
            double h = std::min(dt, burn - elapsed);
            State nx = rk4_step(s, x, h, true);
            if (nx.vg >= vg_end) {
                double lo = 0.0, hi = h;
                for (int i = 0; i < 80; ++i) {
                    const double mid = 0.5 * (lo + hi);
                    if (rk4_step(s, x, mid, true).vg < vg_end) lo = mid; else hi = mid;
                }
                h = 0.5 * (lo + hi);
                nx = rk4_step(s, x, h, true);
                depleted = true;
            }
            x = nx;
            elapsed += h;
            res.p_min = std::min(res.p_min, pressure(s, x));
            res.t_gas_min = std::min(res.t_gas_min, x.temp);
            if (depleted) break;
        }
        rec.duration = elapsed;
        rec.p_end = pressure(s, x);
        rec.t_gas_end = x.temp;
        rec.dv = x.dv - dv0;
        res.fire_time += elapsed;
        res.burns.push_back(rec);
        if (depleted) break;
        double c = 0.0;
        while (c < coast - 1e-12) {
            const double h = std::min(dt, coast - c);
            x = rk4_step(s, x, h, false);
            c += h;
        }
    }
    res.final_state = x;
    return res;
}

// Adaptive Simpson for the independent burn-time quadrature.
inline double simpson_adapt(const std::function<double(double)>& f, double a, double b,
                            double eps, double whole, double fa, double fb, double fm,
                            int depth) {
    const double m = 0.5 * (a + b);
    const double lm = 0.5 * (a + m), rm = 0.5 * (m + b);
    const double flm = f(lm), frm = f(rm);
    const double left = (m - a) / 6.0 * (fa + 4 * flm + fm);
    const double right = (b - m) / 6.0 * (fm + 4 * frm + fb);
    if (depth <= 0 || std::fabs(left + right - whole) <= 15 * eps)
        return left + right + (left + right - whole) / 15.0;
    return simpson_adapt(f, a, m, eps / 2, left, fa, fm, flm, depth - 1) +
           simpson_adapt(f, m, b, eps / 2, right, fm, fb, frm, depth - 1);
}

inline double integrate(const std::function<double(double)>& f, double a, double b,
                        double eps = 1e-10) {
    const double fa = f(a), fb = f(b), fm = f(0.5 * (a + b));
    const double whole = (b - a) / 6.0 * (fa + 4 * fm + fb);
    return simpson_adapt(f, a, b, eps, whole, fa, fb, fm, 50);
}

// Continuous isothermal firing time from Vg0 to depletion:
//     t = integral rho_p dVg / (N mdot(P0 Vg0 / Vg))
inline double isothermal_burn_time(const System& s) {
    auto f = [&](double vg) { return s.rho_p / (s.n_thrusters * s.thr.mdot(s.p0 * s.vg0 / vg)); };
    return integrate(f, s.vg0, s.vg0 + s.v_usable, 1e-9);
}

}  // namespace bd
