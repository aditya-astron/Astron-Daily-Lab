// Astron Built this AP
// Project 033 | Thermal analysis | Minimum-mass pumped-loop radiator panel
// Author: Aditya (aditya-astron)
// Date: 2026-11-12 IST
//
// Verification. Each model is checked against something that does not share
// its code path: a closed form through the incomplete beta function, the
// linearised fin, a finite-difference solver, a brute-force optimum, an RK4
// march of the coolant, and the energy balance.

#include <cmath>
#include <cstdio>

#include "radiator.hpp"

using namespace rad;

static int failures = 0;

static void check(bool ok, const char* what, double got, double want) {
    std::printf("[%s] %-58s got %.10g  ref %.10g\n", ok ? " ok " : "FAIL", what, got, want);
    if (!ok) ++failures;
}

static bool rel_close(double a, double b, double tol) { return std::fabs(a - b) <= tol * std::fabs(b); }

// Incomplete beta B(x; a, b) through the Lentz continued fraction, with the
// symmetry relation used where the fraction converges slowly.
static double betacf(double a, double b, double x) {
    const double tiny = 1e-300;
    double c = 1.0, d = 1.0 - (a + b) * x / (a + 1.0);
    if (std::fabs(d) < tiny) d = tiny;
    d = 1.0 / d;
    double h = d;
    for (int m = 1; m < 1000; ++m) {
        const double m2 = 2.0 * m;
        double aa = m * (b - m) * x / ((a + m2 - 1.0) * (a + m2));
        d = 1.0 + aa * d; if (std::fabs(d) < tiny) d = tiny;
        c = 1.0 + aa / c; if (std::fabs(c) < tiny) c = tiny;
        d = 1.0 / d; h *= d * c;
        aa = -(a + m) * (a + b + m) * x / ((a + m2) * (a + m2 + 1.0));
        d = 1.0 + aa * d; if (std::fabs(d) < tiny) d = tiny;
        c = 1.0 + aa / c; if (std::fabs(c) < tiny) c = tiny;
        d = 1.0 / d;
        const double del = d * c;
        h *= del;
        if (std::fabs(del - 1.0) < 1e-16) break;
    }
    return h;
}

static double beta_complete(double a, double b) {
    return std::exp(std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b));
}

static double beta_incomplete(double x, double a, double b) {
    const double front = std::exp(a * std::log(x) + b * std::log1p(-x));
    if (x < (a + 1.0) / (a + b + 2.0)) return front * betacf(a, b, x) / a;
    return beta_complete(a, b) - front * betacf(b, a, 1.0 - x) / b;
}

// Zero-sink fin, closed form: with theta = T/Tt and w = theta^-5,
//   int_1^thb dtheta / sqrt(theta^5 - 1) = (B(3/10, 1/2) - B(thb^-5; 3/10, 1/2)) / 5
// and J = sqrt(5) Tt^{-3/2} times that integral.
static double J_zero_sink_closed(double Tb, double Tt) {
    const double w0 = std::pow(Tb / Tt, -5.0);
    const double I = (beta_complete(0.3, 0.5) - beta_incomplete(w0, 0.3, 0.5)) / 5.0;
    return std::sqrt(5.0) * std::pow(Tt, -1.5) * I;
}

int main() {
    std::printf("Radiator fin and panel verification\n\n");

    // 1. Shape integral against the incomplete-beta closed form (Ts = 0).
    {
        double worst = 0.0;
        for (double phi : {0.2, 0.5, 0.8, 0.95, 0.999}) {
            const double Tb = 320.0, Tt = phi * Tb;
            const double e = std::fabs(shape_integral(Tb, Tt, 0.0) / J_zero_sink_closed(Tb, Tt) - 1.0);
            worst = std::max(worst, e);
        }
        check(worst < 1e-12, "J(Tb,Tt), Ts=0: max rel. error vs incomplete beta", worst, 0.0);
    }

    // 2. Exact fin round trip: solve_tip inverts L(Tt) and q matches the
    //    closed-form q = sqrt(2 n eps sigma k t (Tb^5 - Tt^5)/5) for Ts = 0.
    {
        const FinSpec f{167.0, 1e-3, 0.85, 2.0, 0.0};
        const double Tb = 350.0, Tt_true = 290.0;
        const double L = std::sqrt(f.k * f.t / (2.0 * f.n_faces * f.eps * SIGMA)) *
                         J_zero_sink_closed(Tb, Tt_true);
        const FinResult r = fin_exact(f, Tb, L);
        check(rel_close(r.Tt, Tt_true, 1e-11), "tip temperature recovered from closed-form L", r.Tt, Tt_true);
        const double q = std::sqrt(2.0 * f.n_faces * f.eps * SIGMA * f.k * f.t *
                                   (std::pow(Tb, 5) - std::pow(Tt_true, 5)) / 5.0);
        check(rel_close(r.q, q, 1e-10), "root heat flow vs closed form (W/m)", r.q, q);
    }

    // 3. Small Tb - Ts: the fin is linear, eta = tanh(mL)/(mL) with
    //    m^2 = 4 n eps sigma Ts^3 / (k t).
    {
        const FinSpec f{167.0, 0.8e-3, 0.8, 1.0, 250.0};
        const double L = 0.12, Tb = 250.01;
        const double m = std::sqrt(4.0 * f.n_faces * f.eps * SIGMA * std::pow(f.Ts, 3) / (f.k * f.t));
        const double eta_lin = std::tanh(m * L) / (m * L);
        const double eta = fin_exact(f, Tb, L).eta;
        check(rel_close(eta, eta_lin, 1e-4), "linear limit: eta vs tanh(mL)/mL", eta, eta_lin);
    }

    // 4. Very high conductivity: isothermal fin, eta -> 1.
    {
        const FinSpec f{1e9, 1e-3, 0.8, 1.0, 220.0};
        const double eta = fin_exact(f, 310.0, 0.1).eta;
        check(std::fabs(eta - 1.0) < 1e-6, "k -> infinity: eta -> 1", eta, 1.0);
    }

    // 5. Finite differences (independent Newton solver) converge to the
    //    exact solution at second order, with a sink, and conserve energy.
    {
        const FinSpec f{167.0, 0.5e-3, 0.8, 1.0, 220.0};
        const double Tb = 318.15, L = 0.10;
        const FinResult ex = fin_exact(f, Tb, L);
        double prev_err = 0.0, order = 0.0, err = 0.0, tip_err = 0.0;
        FinFD fd;
        for (int N : {50, 100, 200, 400}) {
            fd = fin_fd(f, Tb, L, N);
            err = std::fabs(fd.q_radiated - ex.q) / ex.q;
            if (prev_err > 0.0) order = std::log2(prev_err / err);
            prev_err = err;
            tip_err = std::fabs(fd.T.back() - ex.Tt);
        }
        check(err < 2e-6, "FD (N=400) radiated q vs exact, rel. error", err, 0.0);
        check(std::fabs(order - 2.0) < 0.05, "FD convergence order (N 200 -> 400)", order, 2.0);
        check(tip_err < 1e-3, "FD tip temperature vs exact (K)", tip_err, 0.0);
        const double bal = std::fabs(fd.q_conduction - fd.q_radiated) / fd.q_radiated;
        check(bal < 1e-5, "FD energy balance: root conduction vs radiation", bal, 0.0);
    }

    // 6. Minimum-mass single fin. The 1-D optimum over Tt is checked by brute
    //    force: for each thickness the FD solver finds the length that carries
    //    the required q, and the smallest t*L over a fine grid is compared.
    {
        const double q = 50.0, Tb = 300.0, k = 167.0, eps = 0.8, n = 1.0;
        for (double Ts : {0.0, 220.0}) {
            const OptimalFin o = optimal_fin(q, Tb, k, eps, n, Ts);
            double best_A = 1e30, best_t = 0.0;
            for (int i = 0; i <= 300; ++i) {
                const double t = o.t * std::exp(-0.6 + 1.2 * (i + 0.37) / 300.0);  // grid avoids o.t
                const FinSpec f{k, t, eps, n, Ts};
                double lo = 0.0, hi = 20.0 * o.L;
                for (int it = 0; it < 60; ++it) {
                    const double L = 0.5 * (lo + hi);
                    if (fin_fd(f, Tb, L, 400).q_radiated < q) lo = L; else hi = L;
                }
                const double A = t * 0.5 * (lo + hi);
                if (A < best_A) { best_A = A; best_t = t; }
            }
            char msg[96];
            std::snprintf(msg, sizeof msg, "Ts=%3.0f K: optimum fin area vs brute-force FD (m^2)", Ts);
            check(rel_close(o.area, best_A, 2e-5) && o.area <= best_A * (1.0 + 2e-5), msg, o.area, best_A);
            std::snprintf(msg, sizeof msg, "Ts=%3.0f K: optimum thickness vs brute-force grid (m)", Ts);
            check(rel_close(o.t, best_t, 0.02), msg, o.t, best_t);
        }
        // Mass scaling rho * A ~ rho / k, and Tt* independent of q and k.
        const OptimalFin a = optimal_fin(50.0, 300.0, 167.0, 0.8, 1.0, 220.0);
        const OptimalFin b = optimal_fin(80.0, 300.0, 391.0, 0.8, 1.0, 220.0);
        check(std::fabs(a.Tt - b.Tt) < 1e-6, "optimum tip temperature independent of q and k", b.Tt, a.Tt);
        const double ratio = (b.area / a.area) / (std::pow(80.0 / 50.0, 3) * 167.0 / 391.0);
        check(std::fabs(ratio - 1.0) < 1e-9, "optimum area scales as q^3 / k", ratio, 1.0);
    }

    // 7. Panel: tube length by Gauss-Legendre in T against an RK4 march of
    //    dT/dz = -q'(T) / (mdot cp) along the tube, and the energy balance.
    {
        const PanelSpec p{1200.0, 318.15, 288.15, 3600.0, 0.0095, 0.135,
                          2700.0, 167.0, 0.80, 1.0, 220.0};
        const double t = 0.5e-3, L = 0.08;
        const double Lt = tube_length(p, t, L);
        const double mc = p.mdot() * p.cp;
        auto dTdz = [&](double T) { return -q_per_length(p, t, L, T) / mc; };
        double z = 0.0, T = p.Tin;
        const double h = Lt / 2000.0;
        while (true) {
            const double k1 = dTdz(T), k2 = dTdz(T + 0.5 * h * k1);
            const double k3 = dTdz(T + 0.5 * h * k2), k4 = dTdz(T + h * k3);
            const double Tn = T + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
            if (Tn <= p.Tout) {  // final partial step by linear interpolation
                const double frac = (T - p.Tout) / (T - Tn);
                z += frac * h;
                break;
            }
            T = Tn;
            z += h;
        }
        check(rel_close(z, Lt, 1e-6), "tube length: quadrature vs RK4 march (m)", Lt, z);
        // Energy: integral of q' dz along the march must equal Q.
        double Qrad = 0.0;
        T = p.Tin;
        for (double zz = 0.0; zz < Lt - 1e-12;) {
            const double hh = std::min(h, Lt - zz);
            const double k1 = dTdz(T), k2 = dTdz(T + 0.5 * hh * k1);
            const double k3 = dTdz(T + 0.5 * hh * k2), k4 = dTdz(T + hh * k3);
            Qrad += hh / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4) * -mc;
            T += hh / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
            zz += hh;
        }
        check(std::fabs(T - p.Tout) < 1e-4, "outlet temperature after Ltube of march (K)", T, p.Tout);
        check(rel_close(Qrad, p.Q, 1e-6), "heat rejected along the tube (W)", Qrad, p.Q);
        // Real panel area always exceeds the eta = 1 bound.
        const PanelResult r = panel(p, t, L);
        check(r.area > ideal_area(p), "panel area above isothermal bound (m^2)", r.area, ideal_area(p));
    }

    // 8. Panel optimiser: the Nelder-Mead result is stationary and no point
    //    of a surrounding grid is lighter.
    {
        const PanelSpec p{1200.0, 318.15, 288.15, 3600.0, 0.0095, 0.135,
                          2700.0, 167.0, 0.80, 1.0, 220.0};
        const PanelResult o = optimise_panel(p, 1e-3, 0.05, 1e-6);
        double gmin = 1e30;
        for (int i = -10; i <= 10; ++i)
            for (int j = -10; j <= 10; ++j)
                gmin = std::min(gmin, panel(p, o.t * std::exp(0.03 * i), o.L * std::exp(0.03 * j)).mass);
        check(o.mass <= gmin * (1.0 + 1e-9), "optimum mass vs 21x21 surrounding grid (kg)", o.mass, gmin);
        const double d = 1e-3;
        const double gt = (panel(p, o.t * std::exp(d), o.L).mass - panel(p, o.t * std::exp(-d), o.L).mass) / (2 * d);
        const double gL = (panel(p, o.t, o.L * std::exp(d)).mass - panel(p, o.t, o.L * std::exp(-d)).mass) / (2 * d);
        const double gnorm = std::hypot(gt, gL) / o.mass;
        check(gnorm < 1e-4, "optimum stationarity |dM/dln(t,L)| / M", gnorm, 0.0);
    }

    std::printf("\n%s (%d failure%s)\n", failures ? "FAILED" : "ALL PASSED", failures,
                failures == 1 ? "" : "s");
    return failures ? 1 : 0;
}
