// Astron Built this AP
// Project 033 | Thermal analysis | Minimum-mass pumped-loop radiator panel
// Author: Aditya (aditya-astron)
// Date: 2026-11-12 IST
//
// Sizes a body-mounted tube-and-fin radiator that rejects 1.2 kW from a
// single-phase pumped loop (45 degC in, 15 degC out) to a 220 K effective
// sink, and asks how thin and how wide the fins should be for least mass.

#include <cstdio>

#include "radiator.hpp"

using namespace rad;

static PanelSpec baseline_spec(double k = 167.0, double rho = 2700.0) {
    PanelSpec p;
    p.Q = 1200.0;
    p.Tin = 318.15;
    p.Tout = 288.15;
    p.cp = 3600.0;      // 50/50 propylene glycol-water
    p.D = 0.0095;       // 3/8 in tube
    p.m_tube = 0.135;   // tube 0.066 + coolant 0.049 + saddle and bond 0.020 kg/m
    p.rho = rho;
    p.k = k;
    p.eps = 0.80;       // silvered FEP, end of life
    p.n_faces = 1.0;
    p.Ts = 220.0;
    return p;
}

static void print_panel(const char* name, const PanelResult& r) {
    std::printf("  %-26s %6.2f %7.1f %7.1f %7.2f %7.2f %7.2f\n", name, r.t * 1e3,
                (2.0 * r.L + 0.0095) * 1e3, r.Ltube, r.area, r.mass_fins, r.mass);
}

// Least mass with fin gauge >= t_min and planform area <= A_max. Area grows
// with fin length at fixed t, so for each t the active constraint fixes L by
// bisection; golden section then searches t.
static PanelResult min_mass_with_area(const PanelSpec& p, double A_max, double t_min) {
    const PanelResult free_opt = optimise_panel(p, 1e-3, 0.05, t_min);
    if (free_opt.area <= A_max) return free_opt;
    auto at_area = [&](double t) {
        double lo = std::log(1e-3), hi = std::log(0.4);
        for (int it = 0; it < 40; ++it) {
            const double m = 0.5 * (lo + hi);
            if (panel(p, t, std::exp(m)).area > A_max) hi = m; else lo = m;
        }
        return panel(p, t, std::exp(lo));
    };
    const double g = 0.5 * (std::sqrt(5.0) - 1.0);
    double a = std::log(t_min), b = std::log(5e-3);
    double x1 = b - g * (b - a), x2 = a + g * (b - a);
    double f1 = at_area(std::exp(x1)).mass, f2 = at_area(std::exp(x2)).mass;
    for (int it = 0; it < 40; ++it) {
        if (f1 < f2) { b = x2; x2 = x1; f2 = f1; x1 = b - g * (b - a); f1 = at_area(std::exp(x1)).mass; }
        else         { a = x1; x1 = x2; f1 = f2; x2 = a + g * (b - a); f2 = at_area(std::exp(x2)).mass; }
    }
    const PanelResult r = at_area(std::exp(0.5 * (a + b)));
    const PanelResult edge = at_area(t_min);
    return edge.mass < r.mass ? edge : r;
}

int main() {
    const PanelSpec p = baseline_spec();
    const double t_gauge = 0.5e-3;  // thinnest sheet we would bond and handle

    std::printf("Pumped-loop radiator: Q = %.0f W, loop %.1f -> %.1f K, mdot = %.2f g/s\n",
                p.Q, p.Tin, p.Tout, p.mdot() * 1e3);
    std::printf("Al 6061 fins, eps = %.2f, one radiating face, sink %.0f K\n", p.eps, p.Ts);
    std::printf("Isothermal (eta = 1) area bound: %.2f m^2\n\n", ideal_area(p));

    std::printf("Designs                        t mm  pitch mm  tube m   area m2  fins kg  total kg\n");
    print_panel("baseline 1.0 mm, 100 mm", panel(p, 1.0e-3, 0.5 * (0.100 - p.D)));
    const PanelResult free_opt = optimise_panel(p, 1e-3, 0.05, 1e-6);
    print_panel("optimum, no gauge limit", free_opt);
    const PanelResult gauge_opt = optimise_panel(p, 1e-3, 0.05, t_gauge);
    print_panel("optimum, t >= 0.5 mm", gauge_opt);

    std::printf("\nMass vs. planform budget (t >= 0.5 mm)\n");
    std::printf("  A_max m2   t mm  pitch mm  tube m  mass kg   eta_fin at inlet/outlet\n");
    for (double A : {4.6, 4.8, 5.0, 5.5}) {
        const PanelResult r = min_mass_with_area(p, A, t_gauge);
        const FinSpec f = fin_of(p, r.t);
        std::printf("  %7.2f %7.2f %8.1f %7.1f %8.2f        %.3f / %.3f\n", A, r.t * 1e3,
                    (2.0 * r.L + p.D) * 1e3, r.Ltube, r.mass, fin_exact(f, p.Tin, r.L).eta,
                    fin_exact(f, p.Tout, r.L).eta);
    }

    std::printf("\nFin material at t >= 0.5 mm (tube line mass unchanged)\n");
    std::printf("  material        k     rho   k/rho    t mm  pitch mm  area m2  mass kg\n");
    struct Mat { const char* name; double k, rho; };
    for (const Mat& m : {Mat{"Al 1100-O", 222.0, 2710.0}, Mat{"Al 6061-T6", 167.0, 2700.0},
                         Mat{"Mg AZ31B", 96.0, 1770.0}, Mat{"Cu C10100", 391.0, 8940.0}}) {
        const PanelResult r = optimise_panel(baseline_spec(m.k, m.rho), 1e-3, 0.05, t_gauge);
        std::printf("  %-12s %5.0f %7.0f %7.4f %7.2f %8.1f %8.2f %8.2f\n", m.name, m.k, m.rho,
                    m.k / m.rho, r.t * 1e3, (2.0 * r.L + p.D) * 1e3, r.area, r.mass);
    }

    std::printf("\nMinimum-mass single fin (one face, eps 0.80, Tb 300 K, q 50 W/m, Al 6061)\n");
    std::printf("  sink K   Tt/Tb     t mm     L mm    eta\n");
    for (double Ts : {0.0, 150.0, 220.0, 260.0}) {
        const OptimalFin o = optimal_fin(50.0, 300.0, 167.0, 0.80, 1.0, Ts);
        std::printf("  %6.0f %7.4f %8.3f %8.1f %6.3f\n", Ts, o.Tt / 300.0, o.t * 1e3, o.L * 1e3, o.eta);
    }
    return 0;
}
