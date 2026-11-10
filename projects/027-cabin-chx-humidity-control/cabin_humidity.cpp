// Astron Built this AP
// Project 027 | ECLSS | Cabin humidity control with a condensing heat exchanger
// Author: Aditya (aditya-astron)
// Date: 2026-11-10 IST
//
// Steady-state humidity of a crewed module whose CHX loop controls temperature
// only. The bypass valve holds 22 degC; humidity is whatever the coil leaves.
// Questions answered:
//   1. Where does cabin RH / dew point settle for 2..8 crew?
//   2. How does it follow the coolant inlet temperature?
//   3. What is the warmest coolant that still keeps RH <= 60 % with 6 crew,
//      one of them exercising?
#include "chx.hpp"

#include <cstdio>
#include <stdexcept>

using namespace chx;

namespace {

constexpr double T_SET        = 22.0;           // degC cabin setpoint
constexpr double Q_EQUIPMENT  = 1200.0;         // W avionics, lighting, fans into the air
constexpr double Q_CREW_SENS  = 100.0;          // W sensible per crew at rest
constexpr double VAPOR_CREW   = 1.9 / 86400.0;  // kg/s per crew, respiration + perspiration
constexpr double VAPOR_EXER   = 0.35 / 3600.0;  // kg/s extra sweat while exercising
constexpr double Q_EXER_SENS  = 250.0;          // W extra sensible while exercising
constexpr double RH_MAX       = 0.60;

const Loop LOOP{
    {700.0, 1800.0},  // CHX hA air / coolant, W/K
    0.24,             // fan, kg/s dry air (~200 L/s)
    0.075,            // coolant, kg/s
    4.4,              // coolant inlet, degC
};

Loads crew_loads(int crew, int exercising) {
    return {Q_EQUIPMENT + crew * Q_CREW_SENS + exercising * Q_EXER_SENS,
            crew * VAPOR_CREW + exercising * VAPOR_EXER};
}

void print_row(const char* label, const Loads& ld, const CabinState& c) {
    std::printf("%-12s %7.0f %8.2f  %5.1f  %6.2f  %6.3f  %7.2f  %7.2f  %7.0f\n", label,
                ld.q_sensible, ld.m_vapor * 86400.0, 100.0 * c.rh(), c.dew(), c.frac_chx,
                c.hx.out.t, c.hx.t_cool_out, c.hx.q);
}

const char* HEADER =
    "case          Qs[W]  vap[kg/d]  RH[%]  dew[C]  f_chx   Tout[C]  Tcool,o   Qchx[W]\n";

}  // namespace

int main() {
    std::printf("Cabin %.1f degC at %.3f kPa, fan %.3f kg/s, CHX hA %.0f/%.0f W/K, coolant %.3f kg/s\n\n",
                T_SET, P_ATM / 1e3, LOOP.m_fan, LOOP.hx.hA_air, LOOP.hx.hA_cool, LOOP.m_cool);

    std::printf("1) Crew count, coolant inlet %.1f degC\n%s", LOOP.t_cool_in, HEADER);
    for (int n = 2; n <= 8; ++n) {
        char label[16];
        std::snprintf(label, sizeof label, "%d crew", n);
        const Loads ld = crew_loads(n, 0);
        print_row(label, ld, solve_cabin(LOOP, ld, T_SET));
    }
    const Loads peak = crew_loads(6, 1);
    print_row("6 + 1 exer", peak, solve_cabin(LOOP, peak, T_SET));

    std::printf("\n2) Coolant inlet temperature, 4 crew\n%s", HEADER);
    for (double tc = 2.0; tc <= 12.01; tc += 2.0) {
        Loop L = LOOP;
        L.t_cool_in = tc;
        char label[16];
        std::snprintf(label, sizeof label, "Tc,in %4.1f", tc);
        const Loads ld = crew_loads(4, 0);
        try {
            print_row(label, ld, solve_cabin(L, ld, T_SET));
        } catch (const std::runtime_error&) {
            std::printf("%-12s %7.0f %8.2f  setpoint lost: CHX at full flow cannot remove the load\n",
                        label, ld.q_sensible, ld.m_vapor * 86400.0);
        }
    }

    // RH rises monotonically with coolant inlet temperature, so solve RH = RH_MAX.
    // Losing the temperature setpoint counts as a violation.
    auto rh_at = [&](double tc) {
        Loop L = LOOP;
        L.t_cool_in = tc;
        try {
            return solve_cabin(L, peak, T_SET).rh();
        } catch (const std::runtime_error&) {
            return 1.0;
        }
    };
    const double tc_max = find_root([&](double tc) { return rh_at(tc) - RH_MAX; }, 4.0, 14.0, 1e-6);
    Loop L = LOOP;
    L.t_cool_in = tc_max;
    const CabinState c = solve_cabin(L, peak, T_SET);
    std::printf("\n3) Warmest coolant inlet for RH <= %.0f %% with 6 crew, 1 exercising\n%s",
                100.0 * RH_MAX, HEADER);
    print_row("limit", peak, c);
    std::printf("   coolant inlet <= %.2f degC (dew point %.2f degC, CHX carries %.0f %% of fan flow)\n",
                tc_max, c.dew(), 100.0 * c.frac_chx);
    return 0;
}
