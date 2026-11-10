// Astron Built this AP
// Project 027 | ECLSS | Cabin humidity control with a condensing heat exchanger
// Author: Aditya (aditya-astron)
// Date: 2026-11-10 IST
//
// Psychrometrics, a wet-coil condensing heat exchanger (Braun-Klein-Mitchell
// enthalpy effectiveness model plus a segment-by-segment reference model) and
// the steady cabin balance of a CCAA-style loop: one fan, a CHX, and a bypass
// valve that holds cabin temperature. All SI: K or degC as named, Pa, kg/s, W.
#pragma once

#include <cmath>
#include <functional>
#include <stdexcept>

namespace chx {

constexpr double P_ATM   = 101325.0;  // Pa, ISS cabin is held at 14.7 psia
constexpr double CP_DA   = 1006.0;    // J/kg/K, dry air
constexpr double CP_V    = 1860.0;    // J/kg/K, water vapour
constexpr double H_FG0   = 2.501e6;   // J/kg, latent heat at 0 degC
constexpr double CP_W    = 4186.0;    // J/kg/K, liquid coolant (water)
constexpr double EPS_MW  = 0.621945;  // Mw/Ma

// Bracketed root finder: Illinois-modified regula falsi. Keeps the bracket
// like bisection but converges superlinearly, which matters because the cabin
// solve nests four levels of root finding. f(lo) and f(hi) must differ in sign.
inline double find_root(const std::function<double(double)>& f, double lo, double hi,
                     double tol = 1e-12, int max_iter = 200) {
    double flo = f(lo), fhi = f(hi);
    if (flo == 0.0) return lo;
    if (fhi == 0.0) return hi;
    if (flo * fhi > 0.0) throw std::runtime_error("find_root: root not bracketed");
    int side = 0;
    double x = lo, x_prev = hi;
    for (int i = 0; i < max_iter; ++i) {
        x = (lo * fhi - hi * flo) / (fhi - flo);
        if (!(x > lo && x < hi)) x = 0.5 * (lo + hi);   // guard against round-off
        const double fx = f(x);
        if (fx == 0.0 || std::fabs(x - x_prev) < tol) return x;
        x_prev = x;
        if ((fx < 0.0) == (flo < 0.0)) {
            lo = x; flo = fx;
            if (side == -1) fhi *= 0.5;
            side = -1;
        } else {
            hi = x; fhi = fx;
            if (side == +1) flo *= 0.5;
            side = +1;
        }
        if (hi - lo < tol) return x;
    }
    return x;
}

// ---------------------------------------------------------------- psychrometrics

// Saturation pressure over liquid water, Hyland-Wexler (ASHRAE Fundamentals).
// Used below 0 degC as well (supercooled), which is the right branch for a
// condensing surface that never freezes.
inline double p_sat(double t_c) {
    const double T = t_c + 273.15;
    return std::exp(-5.8002206e3 / T + 1.3914993 - 4.8640239e-2 * T
                    + 4.1764768e-5 * T * T - 1.4452093e-8 * T * T * T
                    + 6.5459673 * std::log(T));
}

inline double w_from_pv(double pv, double p = P_ATM) { return EPS_MW * pv / (p - pv); }
inline double pv_from_w(double w, double p = P_ATM) { return p * w / (EPS_MW + w); }
inline double w_sat(double t_c, double p = P_ATM) { return w_from_pv(p_sat(t_c), p); }
inline double w_from_rh(double t_c, double rh, double p = P_ATM) { return w_from_pv(rh * p_sat(t_c), p); }
inline double rh(double t_c, double w, double p = P_ATM) { return pv_from_w(w, p) / p_sat(t_c); }

// Moist-air enthalpy per kg dry air, 0 degC dry air / liquid water reference.
inline double enthalpy(double t_c, double w) { return CP_DA * t_c + w * (H_FG0 + CP_V * t_c); }
inline double temp_from_h(double h, double w) { return (h - w * H_FG0) / (CP_DA + w * CP_V); }
inline double h_sat(double t_c, double p = P_ATM) { return enthalpy(t_c, w_sat(t_c, p)); }
inline double cp_moist(double w) { return CP_DA + w * CP_V; }

inline double dew_point(double w, double p = P_ATM) {
    const double pv = pv_from_w(w, p);
    return find_root([&](double t) { return p_sat(t) - pv; }, -40.0, 80.0, 1e-10);
}
inline double t_from_hsat(double h, double p = P_ATM) {
    return find_root([&](double t) { return h_sat(t, p) - h; }, -40.0, 80.0, 1e-10);
}

// ---------------------------------------------------------------- heat exchanger

struct AirState { double t, w; };   // degC, kg/kg dry air

struct Chx {
    double hA_air;   // W/K, air-side convective conductance (fin efficiency folded in)
    double hA_cool;  // W/K, coolant-side conductance incl. wall
};

struct ChxResult {
    AirState out;
    double q;            // W, total heat removed from the air
    double condensate;   // kg/s
    double t_cool_out;   // degC
    bool wet;
};

inline double eps_counterflow(double ntu, double cr) {
    if (std::fabs(1.0 - cr) < 1e-9) return ntu / (1.0 + ntu);
    const double e = std::exp(-ntu * (1.0 - cr));
    return (1.0 - e) / (1.0 - cr * e);
}

// Dry coil: ordinary counterflow eps-NTU with UA = 1/(1/hA_air + 1/hA_cool).
inline ChxResult chx_dry(const Chx& x, double m_air, AirState in, double m_cool, double t_cool_in) {
    if (x.hA_air <= 0.0 || x.hA_cool <= 0.0) return {in, 0.0, 0.0, t_cool_in, false};
    const double c_a = m_air * cp_moist(in.w), c_w = m_cool * CP_W;
    const double c_min = std::fmin(c_a, c_w), c_max = std::fmax(c_a, c_w);
    const double ua = 1.0 / (1.0 / x.hA_air + 1.0 / x.hA_cool);
    const double q = eps_counterflow(ua / c_min, c_min / c_max) * c_min * (in.t - t_cool_in);
    return {{in.t - q / c_a, in.w}, q, 0.0, t_cool_in + q / c_w, false};
}

// Fully wet coil, Braun-Klein-Mitchell (1989). With Lewis number 1 the air-side
// heat and mass transfer collapse into one enthalpy potential, h_air - h_sat(Ts).
// Linearising h_sat between the coolant inlet and outlet with slope c_s turns
// the coil into a dry exchanger on enthalpy:
//   m*   = m_air c_s / (m_cool cp_w)
//   NTU* = 1 / (m_air (cp_m/hA_air + c_s/hA_cool))
//   Q    = eps*(NTU*, m*) m_air (h_air,in - h_sat(T_cool,in))
// c_s depends on the coolant outlet temperature, so it is iterated.
// The outlet air state comes from the effective surface condition: the air
// approaches a single saturated surface state with NTU_o = hA_air/(m_air cp_m).
inline ChxResult chx_wet(const Chx& x, double m_air, AirState in, double m_cool, double t_cool_in) {
    if (x.hA_air <= 0.0 || x.hA_cool <= 0.0) return {in, 0.0, 0.0, t_cool_in, true};
    const double cpm = cp_moist(in.w);
    const double h_in = enthalpy(in.t, in.w);
    const double hs_ci = h_sat(t_cool_in);
    double t_co = t_cool_in + 1.0, q = 0.0;
    for (int it = 0; it < 100; ++it) {
        const double dt = std::fmax(t_co - t_cool_in, 1e-3);   // secant slope, never 0/0
        const double c_s = (h_sat(t_cool_in + dt) - hs_ci) / dt;
        const double m_star = m_air * c_s / (m_cool * CP_W);
        const double ntu = 1.0 / (m_air * (cpm / x.hA_air + c_s / x.hA_cool));
        // eps* is defined on the smaller capacity. When the coolant's capacity in
        // enthalpy units, m_cool cp_w / c_s, is the smaller one, swap the roles.
        q = (m_star <= 1.0)
            ? eps_counterflow(ntu, m_star) * m_air * (h_in - hs_ci)
            : eps_counterflow(ntu * m_star, 1.0 / m_star) * (m_cool * CP_W / c_s) * (h_in - hs_ci);
        const double t_new = t_cool_in + q / (m_cool * CP_W);
        if (std::fabs(t_new - t_co) < 1e-10) { t_co = t_new; break; }
        t_co = t_new;
    }
    const double h_out = h_in - q / m_air;
    const double ntu_o = x.hA_air / (m_air * cpm);
    const double decay = std::exp(-ntu_o);
    const double h_s_eff = h_in - (h_in - h_out) / (1.0 - decay);
    const double t_s_eff = t_from_hsat(h_s_eff);
    const double w_s_eff = w_sat(t_s_eff);
    AirState out;
    out.w = (in.w > w_s_eff) ? w_s_eff + (in.w - w_s_eff) * decay : in.w;
    out.t = temp_from_h(h_out, out.w);
    return {out, q, m_air * (in.w - out.w), t_co, true};
}

// Dry-surface temperature between local air and coolant temperatures; the
// conductance ratio is the same everywhere because hA is spread uniformly.
inline double t_surface_dry(const Chx& x, double t_air, double t_cool) {
    return t_cool + x.hA_air / (x.hA_air + x.hA_cool) * (t_air - t_cool);
}

// Rating of a counterflow coil that may be dry, fully wet or partly wet.
// Air meets the warm coolant first, so a partly wet coil is dry over the first
// fraction phi of its area and wet over the rest. The boundary sits where the
// dry surface temperature falls to the inlet dew point (w is constant in the
// dry part). For a given phi the two sections are coupled through the coolant
// temperature at the boundary, found by a bracketed root solve; phi is then
// solved from the boundary condition.
inline ChxResult chx_rate(const Chx& x, double m_air, AirState in, double m_cool, double t_cool_in) {
    const double t_dp = dew_point(in.w);
    ChxResult dry = chx_dry(x, m_air, in, m_cool, t_cool_in);
    // Coldest surface: air outlet / coolant inlet end.
    if (t_surface_dry(x, dry.out.t, t_cool_in) >= t_dp) return dry;

    // Warmest wet-model surface: air inlet / coolant outlet end, where
    // hA_cool (Ts - Tc) = hA_air/cp_m (h_in - h_sat(Ts)).
    ChxResult wet = chx_wet(x, m_air, in, m_cool, t_cool_in);
    const double cpm = cp_moist(in.w), h_in = enthalpy(in.t, in.w);
    const double ts_warm = find_root([&](double s) {
        return x.hA_air / cpm * (h_in - h_sat(s)) - x.hA_cool * (s - wet.t_cool_out); },
        wet.t_cool_out, in.t, 1e-10);
    if (ts_warm <= t_dp) return wet;

    struct Split { ChxResult r; double t_boundary_surface; };
    auto split = [&](double phi) {
        const Chx xd{phi * x.hA_air, phi * x.hA_cool};
        const Chx xw{(1.0 - phi) * x.hA_air, (1.0 - phi) * x.hA_cool};
        auto sections = [&](double t_cb, ChxResult& d, ChxResult& w) {
            d = chx_dry(xd, m_air, in, m_cool, t_cb);
            w = chx_wet(xw, m_air, d.out, m_cool, t_cool_in);
        };
        ChxResult d, w;
        const double t_cb = find_root([&](double t) { sections(t, d, w); return w.t_cool_out - t; },
                                   t_cool_in, in.t, 1e-10);
        sections(t_cb, d, w);
        ChxResult r{w.out, d.q + w.q, m_air * (in.w - w.out.w), d.t_cool_out, true};
        return Split{r, t_surface_dry(x, d.out.t, t_cb)};
    };
    const double phi = find_root([&](double p) { return split(p).t_boundary_surface - t_dp; },
                              0.0, 1.0, 1e-12);
    return split(phi).r;
}

// Reference model: N counterflow segments, local surface temperature from the
// series balance hA_cool (Ts - Tc) = air-side flux, Lewis number 1, condensation
// only where w_air > w_sat(Ts). The coolant outlet temperature is found by
// shooting so the coolant enters at t_cool_in. No linearisation of h_sat.
inline ChxResult chx_segmented(const Chx& x, double m_air, AirState in, double m_cool,
                               double t_cool_in, int n = 4000) {
    const double ga = x.hA_air / n, gc = x.hA_cool / n;
    struct March { AirState out; double q, t_cool_end; bool wet; };
    auto march = [&](double t_cool_out) {
        AirState a = in;
        double tc = t_cool_out, q = 0.0;
        bool any_wet = false;
        for (int i = 0; i < n; ++i) {
            const double cpm = cp_moist(a.w);
            const double h_a = enthalpy(a.t, a.w);
            auto flux_at = [&](double ts, bool& wet) {
                wet = a.w > w_sat(ts);
                return wet ? ga / cpm * (h_a - h_sat(ts)) : ga * (a.t - ts);
            };
            bool wet;
            const double ts = find_root([&](double s) { bool w; return flux_at(s, w) - gc * (s - tc); },
                                     std::fmin(tc, a.t) - 1.0, std::fmax(tc, a.t) + 1.0, 1e-9);
            const double dq = flux_at(ts, wet);
            const double dw = wet ? ga / cpm * (a.w - w_sat(ts)) / m_air : 0.0;
            a.w -= dw;
            a.t = temp_from_h(h_a - dq / m_air, a.w);
            tc -= dq / (m_cool * CP_W);   // coolant flows the other way
            q += dq;
            any_wet = any_wet || wet;
        }
        return March{a, q, tc, any_wet};
    };
    const double t_co = find_root([&](double t) { return march(t).t_cool_end - t_cool_in; },
                               t_cool_in, in.t, 1e-8);
    March m = march(t_co);
    return {m.out, m.q, m_air * (in.w - m.out.w), t_co, m.wet};
}

// ---------------------------------------------------------------- cabin balance

struct Loads {
    double q_sensible;   // W, crew sensible + avionics + lighting into the cabin air
    double m_vapor;      // kg/s, crew respiration/perspiration and hygiene water
};

struct Loop {
    Chx hx;
    double m_fan;        // kg/s dry air through the fan
    double m_cool;       // kg/s coolant through the CHX
    double t_cool_in;    // degC
};

struct CabinState {
    double t, w;         // well-mixed cabin
    double frac_chx;     // fraction of fan flow sent through the CHX (rest bypasses)
    ChxResult hx;
    double rh() const { return chx::rh(t, w); }
    double dew() const { return chx::dew_point(w); }
};

// Steady cabin at temperature setpoint t_set. Unknowns: CHX flow fraction f and
// cabin humidity ratio w. Mixing the CHX outlet with the bypass flow and
// balancing the well-mixed cabin gives, per kg dry air,
//   energy:  f m_fan (h_cabin - h_chx,out) = Q_s + m_v h_g(T)
//   water:   f m_fan (w_cabin - w_chx,out) = m_v
// i.e. the CHX must remove the full sensible load plus the vapour enthalpy and
// condense exactly the vapour load. Nested root finding: f for energy at fixed w,
// then w for water.
inline CabinState solve_cabin(const Loop& L, const Loads& ld, double t_set) {
    const double h_g = H_FG0 + CP_V * t_set;
    const double q_need = ld.q_sensible + ld.m_vapor * h_g;
    auto frac_for = [&](double w) {
        auto resid = [&](double f) {
            return chx_rate(L.hx, f * L.m_fan, {t_set, w}, L.m_cool, L.t_cool_in).q - q_need;
        };
        if (resid(1.0) < 0.0) return -1.0;   // CHX cannot hold the setpoint
        return find_root(resid, 1e-6, 1.0, 1e-12);
    };
    // Full CHX flow can only hold the setpoint once the cabin is humid enough
    // for latent removal to make up the difference. Below that, too little
    // water is removed, so report the condensate shortfall as the whole load.
    auto water_resid = [&](double w) {
        const double f = frac_for(w);
        if (f < 0.0) return -ld.m_vapor;
        return chx_rate(L.hx, f * L.m_fan, {t_set, w}, L.m_cool, L.t_cool_in).condensate - ld.m_vapor;
    };
    const double w_lo = w_sat(L.t_cool_in) * 1.0001, w_hi = w_sat(t_set) * 0.999;
    if (water_resid(w_hi) < 0.0)
        throw std::runtime_error("solve_cabin: CHX cannot hold the temperature setpoint");
    const double w = find_root(water_resid, w_lo, w_hi, 1e-13);
    const double f = frac_for(w);
    if (f < 0.0) throw std::runtime_error("solve_cabin: CHX cannot hold the temperature setpoint");
    return {t_set, w, f, chx_rate(L.hx, f * L.m_fan, {t_set, w}, L.m_cool, L.t_cool_in)};
}

}  // namespace chx
