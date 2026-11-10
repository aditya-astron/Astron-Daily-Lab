// Astron Built this AP
// Project 027 | ECLSS | Cabin humidity control with a condensing heat exchanger
// Author: Aditya (aditya-astron)
// Date: 2026-11-10 IST
//
// Verification: psychrometrics against tabulated values, the effectiveness
// models against the segmented reference coil, conservation, and the steady
// cabin solution against a transient simulation with a PI bypass controller.
#include "chx.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>

using namespace chx;

static int failures = 0, checks = 0;

static void check(bool ok, const char* what, double got, double want) {
    ++checks;
    if (!ok) ++failures;
    std::printf("[%s] %-58s got %-14.8g want %.8g\n", ok ? "PASS" : "FAIL", what, got, want);
}
static void close_rel(const char* what, double got, double want, double rtol) {
    check(std::fabs(got - want) <= rtol * std::fabs(want), what, got, want);
}
static void close_abs(const char* what, double got, double want, double atol) {
    check(std::fabs(got - want) <= atol, what, got, want);
}

static void test_psychrometrics() {
    // Saturation pressure over water, IAPWS-95 values (Pa).
    close_rel("p_sat(0.01 C) triple point", p_sat(0.01), 611.657, 1e-3);
    close_rel("p_sat(20 C)", p_sat(20.0), 2339.2, 1e-3);
    close_rel("p_sat(25 C)", p_sat(25.0), 3169.9, 1e-3);
    close_rel("p_sat(40 C)", p_sat(40.0), 7384.9, 1e-3);
    close_rel("p_sat(60 C)", p_sat(60.0), 19946.0, 1e-3);

    // Hand calculation, 22 C 50 % RH at sea level:
    // pv = 0.5 * 2644.7 = 1322.3 Pa, w = 0.621945 pv/(p - pv) = 0.008224
    // h = 1006*22 + w (2.501e6 + 1860*22) = 43.04 kJ/kg (psychrometric chart ~43)
    const double w = w_from_rh(22.0, 0.5);
    close_rel("w at 22 C / 50 %", w, 0.008224, 2e-4);
    close_rel("h at 22 C / 50 % [J/kg]", enthalpy(22.0, w), 43040.0, 1e-3);
    close_abs("dew point at 22 C / 50 % [C] (table 11.1)", dew_point(w), 11.1, 0.05);
    close_abs("dew point round trip", dew_point(w_sat(7.3)), 7.3, 1e-8);
    close_abs("t_from_hsat round trip", t_from_hsat(h_sat(13.7)), 13.7, 1e-8);
    close_abs("rh round trip", rh(22.0, w), 0.5, 1e-12);
}

static void test_dry_coil() {
    // Classical counterflow eps-NTU done by hand for a dry case.
    const Chx hx{700.0, 1800.0};
    const AirState in{22.0, w_from_rh(22.0, 0.3)};   // dew ~3.6 C, coolant 15 C: dry
    const double m_air = 0.12, m_c = 0.075, tci = 15.0;
    const double ca = m_air * (1006.0 + in.w * 1860.0), cw = m_c * 4186.0;
    const double ua = 1.0 / (1.0 / 700.0 + 1.0 / 1800.0), ntu = ua / ca, cr = ca / cw;
    const double eps = (1 - std::exp(-ntu * (1 - cr))) / (1 - cr * std::exp(-ntu * (1 - cr)));
    const ChxResult r = chx_rate(hx, m_air, in, m_c, tci);
    check(!r.wet && r.condensate == 0.0, "dry coil: no condensate", r.condensate, 0.0);
    close_rel("dry coil Q vs hand eps-NTU", r.q, eps * ca * (in.t - tci), 1e-12);
    const ChxResult s = chx_segmented(hx, m_air, in, m_c, tci);
    close_rel("dry coil Q vs segmented reference", r.q, s.q, 1e-3);
}

static void test_wet_coil_vs_reference() {
    // Fully and partly wet coils across flows, coolant temperatures and inlet
    // humidity up to 70 % RH. The effectiveness model linearises h_sat over the
    // wet section; the error grows with latent load and coolant temperature
    // rise, reaching ~1.5 % on Q at 70 % RH. Condensate is compared with a
    // 1e-6 kg/s (0.09 kg/day) floor because near the dry/wet transition it is
    // a small difference of large numbers.
    const Chx hx{700.0, 1800.0};
    double worst_q = 0.0, worst_t = 0.0, worst_m = 0.0;
    int wet_cases = 0;
    for (double m_air : {0.03, 0.08, 0.15, 0.24})
        for (double tci : {3.0, 6.0, 9.0})
            for (double rh_in : {0.35, 0.5, 0.7}) {
                const AirState in{22.0, w_from_rh(22.0, rh_in)};
                const ChxResult a = chx_rate(hx, m_air, in, 0.075, tci);
                const ChxResult b = chx_segmented(hx, m_air, in, 0.075, tci, 3000);
                if (!b.wet) {
                    check(!a.wet || a.condensate < 1e-9, "dry reference case also dry", a.condensate, 0.0);
                    continue;
                }
                ++wet_cases;
                worst_q = std::max(worst_q, std::fabs(a.q - b.q) / b.q);
                worst_t = std::max(worst_t, std::fabs(a.out.t - b.out.t));
                worst_m = std::max(worst_m, std::fabs(a.condensate - b.condensate) / (0.01 * b.condensate + 1e-6));
            }
    check(wet_cases >= 25, "wet or partly wet reference cases", wet_cases, 25);
    check(worst_q < 0.02, "worst |dQ|/Q vs segmented", worst_q, 0.02);
    check(worst_t < 0.3, "worst |dT_out| [K] vs segmented", worst_t, 0.3);
    check(worst_m <= 1.0, "worst condensate error / (1% + 1e-6 kg/s)", worst_m, 1.0);

    // At the actual cabin operating point (4 crew, 39 % RH) the model is much tighter.
    const Loop L{{700.0, 1800.0}, 0.24, 0.075, 4.4};
    const CabinState c = solve_cabin(L, {1600.0, 4 * 1.9 / 86400.0}, 22.0);
    const ChxResult ref = chx_segmented(L.hx, c.frac_chx * L.m_fan, {c.t, c.w}, L.m_cool, L.t_cool_in);
    close_rel("operating point: Q vs segmented", c.hx.q, ref.q, 2e-3);
    close_rel("operating point: condensate vs segmented", c.hx.condensate, ref.condensate, 5e-3);
}

static void test_conservation() {
    const Chx hx{700.0, 1800.0};
    const AirState in{22.0, w_from_rh(22.0, 0.55)};
    for (double m_air : {0.05, 0.2}) {
        const ChxResult r = chx_rate(hx, m_air, in, 0.075, 5.0);
        close_rel("energy: Q = m_air (h_in - h_out)", r.q,
                  m_air * (enthalpy(in.t, in.w) - enthalpy(r.out.t, r.out.w)), 1e-9);
        close_rel("energy: Q = m_cool cp dT_cool", r.q, 0.075 * CP_W * (r.t_cool_out - 5.0), 1e-9);
        check(rh(r.out.t, r.out.w) <= 1.0 + 1e-9, "outlet air not supersaturated", rh(r.out.t, r.out.w), 1.0);
    }
}

// Transient cabin: well-mixed air mass M, bypass fraction set by a PI loop on
// cabin temperature (as the CCAA temperature control valve does). Integrate
// to steady state and compare with the algebraic solution.
static void test_cabin_vs_transient() {
    const Loop L{{700.0, 1800.0}, 0.24, 0.075, 4.4};
    const Loads ld{1600.0, 4 * 1.9 / 86400.0};
    const double t_set = 22.0;
    const CabinState ss = solve_cabin(L, ld, t_set);

    const double h_g = H_FG0 + CP_V * t_set;
    const double cond = ss.hx.condensate;
    close_rel("steady: condensate = vapour load", cond, ld.m_vapor, 1e-6);
    close_rel("steady: CHX duty = Qs + m_v h_g", ss.hx.q, ld.q_sensible + ld.m_vapor * h_g, 1e-6);
    check(ss.frac_chx > 0.0 && ss.frac_chx < 1.0, "steady: bypass fraction inside (0,1)", ss.frac_chx, 0.5);

    const double M = 100.0 * 1.19;           // 100 m^3 module
    const double kp = 0.15, ki = kp / 300.0; // 1/K, 1/(K s)
    double t = 25.0, w = w_from_rh(25.0, 0.7), integ = 0.4 / ki, f = 0.4;
    const double dt = 2.0;
    for (int step = 0; step < 48 * 3600 / 2; ++step) {
        const double err = t - t_set;
        f = std::clamp(kp * err + ki * integ, 1e-4, 1.0);
        if (f > 1e-4 && f < 1.0) integ += err * dt;   // conditional integration anti-windup
        const ChxResult r = chx_rate(L.hx, f * L.m_fan, {t, w}, L.m_cool, L.t_cool_in);
        const double mdot = f * L.m_fan;
        const double h = enthalpy(t, w);
        // Mixing the CHX outlet back in removes mdot (h - h_out) and mdot (w - w_out).
        const double dh = (ld.q_sensible + ld.m_vapor * (H_FG0 + CP_V * t)
                           - mdot * (h - enthalpy(r.out.t, r.out.w))) / M;
        const double dw = (ld.m_vapor - mdot * (w - r.out.w)) / M;
        w += dw * dt;
        t = temp_from_h(h + dh * dt, w);
    }
    close_abs("transient: cabin temperature settles on setpoint [C]", t, t_set, 1e-4);
    close_rel("transient: w matches steady solution", w, ss.w, 1e-4);
    close_rel("transient: bypass fraction matches steady solution", f, ss.frac_chx, 1e-3);
}

static void test_trends() {
    const Loop L{{700.0, 1800.0}, 0.24, 0.075, 4.4};
    double prev = 0.0;
    bool mono = true;
    for (int n = 1; n <= 8; ++n) {
        const double r = solve_cabin(L, {1200.0 + 100.0 * n, n * 1.9 / 86400.0}, 22.0).rh();
        mono = mono && r > prev;
        prev = r;
    }
    check(mono, "cabin RH rises with crew count", prev, 0.0);
    // The cabin dew point can never fall below the coolant inlet temperature.
    const CabinState c = solve_cabin(L, {1400.0, 2 * 1.9 / 86400.0}, 22.0);
    check(c.dew() > L.t_cool_in, "cabin dew point above coolant inlet", c.dew(), L.t_cool_in);
}

int main() {
    test_psychrometrics();
    test_dry_coil();
    test_wet_coil_vs_reference();
    test_conservation();
    test_cabin_vs_transient();
    test_trends();
    std::printf("\n%d/%d checks passed\n", checks - failures, checks);
    return failures ? 1 : 0;
}
