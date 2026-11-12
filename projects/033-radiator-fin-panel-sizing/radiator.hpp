// Astron Built this AP
// Project 033 | Thermal analysis | Minimum-mass pumped-loop radiator panel
// Author: Aditya (aditya-astron)
// Date: 2026-11-12 IST
//
// Models for a tube-and-fin space radiator.
//
// Fin: straight rectangular plate of thickness t and length L, root held at
// the tube wall temperature Tb, adiabatic tip, radiating from n faces with
// emissivity eps to an effective sink Ts. Per unit depth:
//
//     k t T'' = n eps sigma (T^4 - Ts^4),   T(0) = Tb,  T'(L) = 0
//
// Multiplying by T' and integrating once from the tip gives the exact first
// integral, with a = n eps sigma / (k t) and Tt the tip temperature:
//
//     T'^2 = 2 a H(T, Tt),   H(T, Tt) = (T^5 - Tt^5)/5 - Ts^4 (T - Tt)
//
// so the root heat flow is q = k t sqrt(2 a H(Tb, Tt)) and the length is a
// single quadrature, L = J(Tb, Tt) / sqrt(2a) with J = int dT / sqrt(H).
// Writing H = (T - Tt) P(T, Tt) and T = Tt + u^2 removes the square-root
// singularity at the tip, so J is a smooth integral in u.
//
// Panel: fluid flows along each tube; one tube carries two fins (half the tube
// pitch each side) and its own footprint D also radiates. Wall temperature is
// taken equal to the local fluid temperature (fluid film and bond resistances
// neglected). The tube length needed to cool the loop from Tin to Tout is
//
//     Ltube = int_{Tout}^{Tin} mdot cp / q'(T) dT
//
// where q'(T) is the heat rejected per metre of tube.

#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

namespace rad {

constexpr double SIGMA = 5.670374419e-8;  // W m^-2 K^-4 (CODATA 2018, exact SI)

struct FinSpec {
    double k;        // thermal conductivity, W/(m K)
    double t;        // thickness, m
    double eps;      // emissivity of radiating faces
    double n_faces;  // 1 = body-mounted panel, 2 = free-standing fin
    double Ts;       // effective sink temperature, K
};

// P(T, Tt) = H(T, Tt) / (T - Tt), written without cancellation.
inline double P(double T, double Tt, double Ts) {
    const double s = T * T * T * T + T * T * T * Tt + T * T * Tt * Tt + T * Tt * Tt * Tt +
                     Tt * Tt * Tt * Tt;
    return s / 5.0 - Ts * Ts * Ts * Ts;
}

inline double H(double T, double Tt, double Ts) { return (T - Tt) * P(T, Tt, Ts); }

// J(Tb, Tt) = int_{Tt}^{Tb} dT / sqrt(H(T, Tt)). With T = Tt + u^2 the
// integrand 2 / sqrt(P) is bounded, but when Tt approaches Ts it peaks in a
// layer of width ~ sqrt(Tt - Ts) at u = 0. u = c sinh(s) with c = sqrt(Tt - Ts)
// stretches that layer, leaving a smooth integrand for Gauss-Legendre.
inline double gauss24(const std::function<double(double)>& f, double a, double b);

inline double shape_integral(double Tb, double Tt, double Ts) {
    if (Tb <= Tt) return 0.0;
    const double c = std::sqrt(Tt - Ts), smax = std::asinh(std::sqrt(Tb - Tt) / c);
    auto f = [&](double s) {
        const double u = c * std::sinh(s);
        return 2.0 * c * std::cosh(s) / std::sqrt(P(Tt + u * u, Tt, Ts));
    };
    const int panels = 6;
    double J = 0.0;
    for (int i = 0; i < panels; ++i) J += gauss24(f, smax * i / panels, smax * (i + 1) / panels);
    return J;
}

inline double fin_a(const FinSpec& f) { return f.n_faces * f.eps * SIGMA / (f.k * f.t); }

inline double fin_length_from_tip(const FinSpec& f, double Tb, double Tt) {
    return shape_integral(Tb, Tt, f.Ts) / std::sqrt(2.0 * fin_a(f));
}

inline double fin_heat_from_tip(const FinSpec& f, double Tb, double Tt) {
    return f.k * f.t * std::sqrt(2.0 * fin_a(f) * H(Tb, Tt, f.Ts));
}

// Tip temperature of a fin of length L. L(Tt) falls monotonically from
// infinity (Tt -> Ts) to zero (Tt -> Tb). With x = ln((Tt - Ts)/(Tb - Ts)),
// ln L(x) is smooth and monotone, so the Illinois variant of
// regula falsi brackets the root and converges superlinearly.
inline double solve_tip(const FinSpec& f, double Tb, double L) {
    if (L <= 0.0) return Tb;
    const double span = Tb - f.Ts;
    auto Tt_of = [&](double x) { return f.Ts + span * std::exp(x); };
    auto g = [&](double x) { return std::log(fin_length_from_tip(f, Tb, Tt_of(x)) / L); };
    // Lower bracket: Tt - Ts no smaller than 1e-9 Ts, where P(Tt, Tt) = Tt^4 - Ts^4
    // still has ~8 good digits. A fin longer than that is infinitely long in
    // practice (the linear analogue is mL > 20).
    const double e_min = std::max(span * std::exp(-25.0), 1e-9 * f.Ts);
    double a = std::log(e_min / span), b = std::log1p(-1e-9);
    if (a >= b) return Tb;
    double ga = g(a), gb = g(b);
    if (ga <= 0.0) return Tt_of(a);  // longer than any fin we can resolve
    if (gb >= 0.0) return Tt_of(b);
    int side = 0;
    for (int it = 0; it < 200; ++it) {
        const double x = (a * gb - b * ga) / (gb - ga);
        const double gx = g(x);
        if (gx > 0.0) { a = x; ga = gx; if (side == -1) gb *= 0.5; side = -1; }
        else          { b = x; gb = gx; if (side == +1) ga *= 0.5; side = +1; }
        if (std::fabs(gx) < 1e-14 || b - a < 1e-14) return Tt_of(x);
    }
    return Tt_of(0.5 * (a + b));
}

struct FinResult {
    double q;    // root heat flow per unit depth, W/m
    double Tt;   // tip temperature, K
    double eta;  // fin efficiency: q / (n eps sigma L (Tb^4 - Ts^4))
};

inline FinResult fin_exact(const FinSpec& f, double Tb, double L) {
    const double Tt = solve_tip(f, Tb, L);
    const double q = fin_heat_from_tip(f, Tb, Tt);
    const double q_ideal = f.n_faces * f.eps * SIGMA * L * (std::pow(Tb, 4) - std::pow(f.Ts, 4));
    return {q, Tt, q / q_ideal};
}

// Independent finite-difference solution: N intervals, Newton on the
// nonlinear tridiagonal system, ghost node for the adiabatic tip.
struct FinFD {
    std::vector<double> T;
    double q_conduction;  // -k t T'(0), second-order one-sided difference
    double q_radiated;    // trapezoidal integral of the surface loss
};

inline FinFD fin_fd(const FinSpec& f, double Tb, double L, int N) {
    const double h = L / N, c = f.k * f.t / (h * h);
    const double r = f.n_faces * f.eps * SIGMA, Ts4 = std::pow(f.Ts, 4);
    std::vector<double> T(N + 1, Tb), F(N + 1), lo(N + 1), di(N + 1), up(N + 1), dx(N + 1);
    for (int it = 0; it < 50; ++it) {
        for (int i = 1; i <= N; ++i) {
            const double Tl = T[i - 1], Tr = (i < N) ? T[i + 1] : T[N - 1];
            F[i] = c * (Tl - 2.0 * T[i] + Tr) - r * (std::pow(T[i], 4) - Ts4);
            lo[i] = (i < N) ? c : 2.0 * c;
            di[i] = -2.0 * c - 4.0 * r * std::pow(T[i], 3);
            up[i] = c;
        }
        // Thomas algorithm on J dx = -F for unknowns 1..N (T[0] is fixed).
        std::vector<double> cp(N + 1), dp(N + 1);
        cp[1] = up[1] / di[1];
        dp[1] = -F[1] / di[1];
        for (int i = 2; i <= N; ++i) {
            const double m = di[i] - lo[i] * cp[i - 1];
            cp[i] = up[i] / m;
            dp[i] = (-F[i] - lo[i] * dp[i - 1]) / m;
        }
        dx[N] = dp[N];
        for (int i = N - 1; i >= 1; --i) dx[i] = dp[i] - cp[i] * dx[i + 1];
        double step = 0.0;
        for (int i = 1; i <= N; ++i) {
            T[i] += dx[i];
            step = std::max(step, std::fabs(dx[i]));
        }
        if (step < 1e-11 * Tb) break;
    }
    FinFD out;
    out.T = T;
    out.q_conduction = f.k * f.t * (3.0 * T[0] - 4.0 * T[1] + T[2]) / (2.0 * h);
    double s = 0.0;
    for (int i = 0; i <= N; ++i) {
        const double w = (i == 0 || i == N) ? 0.5 : 1.0;
        s += w * r * (std::pow(T[i], 4) - Ts4);
    }
    out.q_radiated = s * h;
    return out;
}

// ---------------------------------------------------------------- panel ---

struct PanelSpec {
    double Q;             // heat to reject, W
    double Tin, Tout;     // loop temperatures at radiator inlet and outlet, K
    double cp;            // coolant specific heat, J/(kg K)
    double D;             // tube footprint width that radiates, m
    double m_tube;        // tube + coolant + bond mass per metre of tube, kg/m
    double rho;           // fin material density, kg/m^3
    double k;             // fin conductivity, W/(m K)
    double eps;           // emissivity
    double n_faces;       // radiating faces of the panel
    double Ts;            // effective sink temperature, K
    double mdot() const { return Q / (cp * (Tin - Tout)); }
};

inline FinSpec fin_of(const PanelSpec& p, double t) { return {p.k, t, p.eps, p.n_faces, p.Ts}; }

// Heat rejected per metre of tube at wall temperature T: two fins of length
// L plus the radiating tube footprint.
inline double q_per_length(const PanelSpec& p, double t, double L, double T) {
    const FinSpec f = fin_of(p, t);
    const double q_fins = 2.0 * fin_exact(f, T, L).q;
    const double q_tube = p.n_faces * p.eps * SIGMA * p.D * (std::pow(T, 4) - std::pow(p.Ts, 4));
    return q_fins + q_tube;
}

// Gauss-Legendre 24-point rule on [a, b].
inline double gauss24(const std::function<double(double)>& f, double a, double b) {
    static const double x[12] = {0.0640568928626056, 0.1911188674736163, 0.3150426796961634,
                                 0.4337935076260451, 0.5454214713888395, 0.6480936519369756,
                                 0.7401241915785544, 0.8200019859739029, 0.8864155270044011,
                                 0.9382745520027328, 0.9747285559713095, 0.9951872199970214};
    static const double w[12] = {0.1279381953467522, 0.1258374563468283, 0.1216704729278034,
                                 0.1155056680537256, 0.1074442701159656, 0.0976186521041139,
                                 0.0861901615319533, 0.0733464814110803, 0.0592985849154368,
                                 0.0442774388174198, 0.0285313886289337, 0.0123412297999872};
    const double c = 0.5 * (a + b), h = 0.5 * (b - a);
    double s = 0.0;
    for (int i = 0; i < 12; ++i) s += w[i] * (f(c - h * x[i]) + f(c + h * x[i]));
    return s * h;
}

inline double tube_length(const PanelSpec& p, double t, double L) {
    const double mc = p.mdot() * p.cp;
    return gauss24([&](double T) { return mc / q_per_length(p, t, L, T); }, p.Tout, p.Tin);
}

struct PanelResult {
    double t, L;        // fin thickness and fin length (half the free gap), m
    double Ltube;       // total tube length, m
    double area;        // planform area, m^2
    double mass;        // tubes + fins, kg
    double mass_fins;   // fins only, kg
};

inline PanelResult panel(const PanelSpec& p, double t, double L) {
    const double Lt = tube_length(p, t, L);
    const double mf = Lt * 2.0 * p.rho * t * L;
    return {t, L, Lt, Lt * (p.D + 2.0 * L), Lt * p.m_tube + mf, mf};
}

// Lower bound on area: whole panel at the local fluid temperature (eta = 1).
inline double ideal_area(const PanelSpec& p) {
    const double mc = p.mdot() * p.cp, r = p.n_faces * p.eps * SIGMA, Ts4 = std::pow(p.Ts, 4);
    return gauss24([&](double T) { return mc / (r * (std::pow(T, 4) - Ts4)); }, p.Tout, p.Tin);
}

// Nelder-Mead on (ln t, ln L), with t held at or above a minimum gauge.
inline PanelResult optimise_panel(const PanelSpec& p, double t0, double L0, double t_min) {
    auto obj = [&](const double* x) {
        const double t = std::max(std::exp(x[0]), t_min);
        return panel(p, t, std::exp(x[1])).mass;
    };
    double s[3][2] = {{std::log(t0), std::log(L0)},
                      {std::log(t0) + 0.3, std::log(L0)},
                      {std::log(t0), std::log(L0) + 0.3}};
    double fs[3];
    for (int i = 0; i < 3; ++i) fs[i] = obj(s[i]);
    for (int it = 0; it < 400; ++it) {
        int o[3] = {0, 1, 2};
        std::sort(o, o + 3, [&](int a, int b) { return fs[a] < fs[b]; });
        const int b = o[0], m = o[1], w = o[2];
        if (std::fabs(fs[w] - fs[b]) < 1e-12 * fs[b] &&
            std::hypot(s[w][0] - s[b][0], s[w][1] - s[b][1]) < 1e-9)
            break;
        double c[2], xr[2], xe[2], xc[2];
        for (int d = 0; d < 2; ++d) {
            c[d] = 0.5 * (s[b][d] + s[m][d]);
            xr[d] = c[d] + (c[d] - s[w][d]);
        }
        const double fr = obj(xr);
        if (fr < fs[b]) {
            for (int d = 0; d < 2; ++d) xe[d] = c[d] + 2.0 * (c[d] - s[w][d]);
            const double fe = obj(xe);
            if (fe < fr) { s[w][0] = xe[0]; s[w][1] = xe[1]; fs[w] = fe; }
            else         { s[w][0] = xr[0]; s[w][1] = xr[1]; fs[w] = fr; }
        } else if (fr < fs[m]) {
            s[w][0] = xr[0]; s[w][1] = xr[1]; fs[w] = fr;
        } else {
            for (int d = 0; d < 2; ++d) xc[d] = c[d] + 0.5 * (s[w][d] - c[d]);
            const double fc = obj(xc);
            if (fc < fs[w]) { s[w][0] = xc[0]; s[w][1] = xc[1]; fs[w] = fc; }
            else {
                for (int i : {m, w}) {
                    for (int d = 0; d < 2; ++d) s[i][d] = s[b][d] + 0.5 * (s[i][d] - s[b][d]);
                    fs[i] = obj(s[i]);
                }
            }
        }
    }
    int b = 0;
    for (int i = 1; i < 3; ++i) if (fs[i] < fs[b]) b = i;
    return panel(p, std::max(std::exp(s[b][0]), t_min), std::exp(s[b][1]));
}

// ------------------------------------------------- minimum-mass single fin --
//
// For a required q, eliminating t and L between q = sqrt(2 n eps sigma k t H)
// and L = sqrt(k t / (2 n eps sigma)) J gives the profile area
//
//     A = t L = q^3 J(Tt) / ((2 n eps sigma)^2 k H(Tb, Tt)^{3/2})
//
// so the optimum tip temperature minimises J / H^{3/2} alone. It does not
// depend on q or k, and the minimum fin mass rho A scales exactly as rho / k.

inline double fin_shape_cost(double Tb, double Tt, double Ts) {
    return shape_integral(Tb, Tt, Ts) / std::pow(H(Tb, Tt, Ts), 1.5);
}

struct OptimalFin { double Tt, t, L, area, eta; };

inline OptimalFin optimal_fin(double q, double Tb, double k, double eps, double n, double Ts) {
    // Golden section on Tt in (Ts, Tb); the cost is unimodal there.
    const double g = 0.5 * (std::sqrt(5.0) - 1.0);
    double a = Ts + 1e-6 * (Tb - Ts), b = Tb - 1e-6 * (Tb - Ts);
    double x1 = b - g * (b - a), x2 = a + g * (b - a);
    double f1 = fin_shape_cost(Tb, x1, Ts), f2 = fin_shape_cost(Tb, x2, Ts);
    for (int it = 0; it < 120 && b - a > 1e-12 * Tb; ++it) {
        if (f1 < f2) { b = x2; x2 = x1; f2 = f1; x1 = b - g * (b - a); f1 = fin_shape_cost(Tb, x1, Ts); }
        else         { a = x1; x1 = x2; f1 = f2; x2 = a + g * (b - a); f2 = fin_shape_cost(Tb, x2, Ts); }
    }
    const double Tt = 0.5 * (a + b), r2 = 2.0 * n * eps * SIGMA, Hb = H(Tb, Tt, Ts);
    const double t = q * q / (r2 * k * Hb);
    const double L = std::sqrt(k * t / r2) * shape_integral(Tb, Tt, Ts);
    const double eta = q / (n * eps * SIGMA * L * (std::pow(Tb, 4) - std::pow(Ts, 4)));
    return {Tt, t, L, t * L, eta};
}

}  // namespace rad
