// Astron Built this AP
// Project 037 | ECLSS | Water recovery sorbent bed breakthrough and LUB scale-up
// Author: Aditya (aditya-astron)
// Date: 2026-11-14 IST
//
// Fixed-bed adsorption of a dissolved organic contaminant in a water processor
// polishing bed. Plug flow, linear driving force (LDF) uptake, Langmuir isotherm:
//
//   eps dc/dt + eps u dc/dz + rho_b dq/dt = 0
//   dq/dt = k (q*(c) - q),      q*(c) = qm b c / (1 + b c)
//
// With the retarded time tau = t - z/u the pair becomes two first-order ODEs,
// one along each characteristic direction:
//
//   dc/dz |tau = -(rho_b k / (eps u)) (q*(c) - q)
//   dq/dtau|z  =  k (q*(c) - q)
//
// with c(0, tau) = c0 and q(z, 0) = 0 (bed clean when the front first arrives).
// Both directions are discretised with the trapezoidal rule, which is
// second-order and unconditionally stable here. Each node is one monotone
// scalar equation in c, solved by Newton.

#pragma once

#include <cmath>
#include <stdexcept>
#include <vector>

namespace bed {

struct Isotherm {
    double qm;   // monolayer capacity [mg/g]
    double b;    // Langmuir affinity [cm^3/mg]; b -> 0 with qm*b fixed is linear
    bool linear; // use q* = K c with K = qm*b
    double q(double c) const { return linear ? qm * b * c : qm * b * c / (1.0 + b * c); }
    double dq(double c) const {
        if (linear) return qm * b;
        const double d = 1.0 + b * c;
        return qm * b / (d * d);
    }
};

struct Bed {
    double length;  // [cm]
    double eps;     // bed void fraction [-]
    double rho_b;   // bulk sorbent density [g/cm^3]
    double u;       // interstitial velocity [cm/h]
    double k;       // LDF rate constant [1/h]
    double c0;      // feed concentration [mg/cm^3]
    Isotherm iso;
};

struct Breakthrough {
    std::vector<double> tau;  // retarded time at the outlet [h]
    std::vector<double> x;    // c_out / c0
};

// Stoichiometric (equilibrium) breakthrough in retarded time: the time at which
// the feed delivered equals the bed's equilibrium capacity, minus the plug-flow
// transit L/u. Exact for any kinetics by an overall mass balance.
inline double tau_stoich(const Bed& B) {
    return B.rho_b * B.iso.q(B.c0) * B.length / (B.eps * B.u * B.c0);
}

inline Breakthrough simulate(const Bed& B, int nz, double tau_end, int ntau) {
    if (nz < 1 || ntau < 1) throw std::invalid_argument("grid");
    const double dz = B.length / nz;
    const double dt = tau_end / ntau;
    const double a = 0.5 * B.k * dt;                               // tau-direction weight
    const double bz = 0.5 * dz * B.rho_b * B.k / (B.eps * B.u);    // z-direction weight
    const Isotherm& f = B.iso;
    // Clean-bed transfer units per z cell, 2*bz*q*'(0). The trapezoid does not damp
    // stiff modes, so beyond ~2 the leading edge of the front rings and leaks.
    if (2.0 * bz * f.dq(0.0) > 2.0) throw std::invalid_argument("z grid too coarse for the clean-bed uptake rate");

    // Solve the node equation for c given its two upstream neighbours.
    // cL,qL: same z, previous tau. cB,rB: previous z, same tau (rB = q*-q there).
    auto node = [&](double cL, double qL, double cB, double rB, double& c, double& q) {
        const double A = qL + a * (f.q(cL) - qL);
        // g(c) = c - cB + bz (rB + (q*(c) - A)/(1+a)), g' > 0 and g convex-up in c
        c = cB;
        for (int it = 0; it < 60; ++it) {
            const double g = c - cB + bz * (rB + (f.q(c) - A) / (1.0 + a));
            const double gp = 1.0 + bz * f.dq(c) / (1.0 + a);
            const double step = g / gp;
            c -= step;
            if (c < 0.0) c = 0.0;
            if (std::fabs(step) <= 1e-15 * (B.c0 + c)) break;
        }
        q = (A + a * f.q(c)) / (1.0 + a);
    };

    std::vector<double> c(nz + 1), q(nz + 1, 0.0);
    // tau = 0: clean bed, c decays along z under first-front uptake.
    c[0] = B.c0;
    // With q = 0 the z-trapezoid is cn = c_i - bz (q*(c_i) + q*(cn)).
    for (int i = 0; i < nz; ++i) {
        double cn = c[i];
        for (int it = 0; it < 60; ++it) {
            const double g = cn - c[i] + bz * (f.q(c[i]) + f.q(cn));
            const double step = g / (1.0 + bz * f.dq(cn));
            cn -= step;
            if (cn < 0.0) cn = 0.0;
            if (std::fabs(step) <= 1e-15 * (B.c0 + cn)) break;
        }
        c[i + 1] = cn;
    }

    Breakthrough out;
    out.tau.reserve(ntau + 1);
    out.x.reserve(ntau + 1);
    out.tau.push_back(0.0);
    out.x.push_back(c[nz] / B.c0);

    for (int j = 0; j < ntau; ++j) {
        // Inlet: c = c0, q relaxes towards q*(c0).
        const double qs0 = f.q(B.c0);
        const double qin = (q[0] * (1.0 - a) + 2.0 * a * qs0) / (1.0 + a);
        double cB = B.c0, qB = qin;
        c[0] = B.c0;
        q[0] = qin;
        for (int i = 1; i <= nz; ++i) {
            double cn, qn;
            node(c[i], q[i], cB, f.q(cB) - qB, cn, qn);
            c[i] = cn;
            q[i] = qn;
            cB = cn;
            qB = qn;
        }
        out.tau.push_back((j + 1) * dt);
        out.x.push_back(c[nz] / B.c0);
    }
    return out;
}

// Retarded time at which the outlet first reaches fraction x (linear interpolation).
inline double time_at(const Breakthrough& bt, double x) {
    for (size_t j = 1; j < bt.x.size(); ++j)
        if (bt.x[j] >= x && bt.x[j - 1] < x) {
            const double s = (x - bt.x[j - 1]) / (bt.x[j] - bt.x[j - 1]);
            return bt.tau[j - 1] + s * (bt.tau[j] - bt.tau[j - 1]);
        }
    return NAN;
}

// First moment of the step response, int (1 - x) dtau (trapezoidal). The run
// must extend until the outlet has saturated.
inline double first_moment(const Breakthrough& bt) {
    double s = 0.0;
    for (size_t j = 1; j < bt.x.size(); ++j)
        s += 0.5 * (bt.tau[j] - bt.tau[j - 1]) * ((1.0 - bt.x[j]) + (1.0 - bt.x[j - 1]));
    return s;
}

// Constant-pattern asymptote for Langmuir + solid LDF (no dispersion).
// In a fully developed wave q/q0 = c/c0 = x, and dx/dt = k (1-R) x(1-x) / (R + (1-R) x),
// R = 1/(1 + b c0). Integrating and fixing the constant by the mass balance
// int (x - H(t - t_st)) dt = 0 gives
//   tau(x) = tau_st - 1/k + [R ln x - ln(1 - x)] / (k (1 - R)).
inline double constant_pattern_tau(const Bed& B, double x) {
    const double R = 1.0 / (1.0 + B.iso.b * B.c0);
    return tau_stoich(B) - 1.0 / B.k + (R * std::log(x) - std::log1p(-x)) / (B.k * (1.0 - R));
}

// Anzelius/Thomas solution for a linear isotherm with LDF:
//   x = J(xi, eta) = 1 - int_0^xi exp(-eta - s) I0(2 sqrt(eta s)) ds
//   xi = rho_b k K L / (eps u) (transfer units), eta = k tau.
// The integrand is written as exp(-(sqrt(eta)-sqrt(s))^2) * [exp(-y) I0(y)],
// y = 2 sqrt(eta s), so nothing overflows; composite Gauss-Legendre in sqrt(s)
// removes the endpoint behaviour.
inline double bessel_i0_scaled(double y) {
    // exp(-y) I0(y). std::cyl_bessel_i is fine below ~700; use the asymptotic
    // series above 50 where it is accurate to double precision.
    if (y < 50.0) return std::exp(-y) * std::cyl_bessel_i(0.0, y);
    const double t = 1.0 / (8.0 * y);
    double term = 1.0, sum = 1.0;
    for (int n = 1; n < 30; ++n) {
        term *= (2.0 * n - 1) * (2.0 * n - 1) * t / n;
        sum += term;
        if (term < 1e-17 * sum) break;
    }
    return sum / std::sqrt(2.0 * M_PI * y);
}

inline double thomas_J(double xi, double eta) {
    // s = w^2, ds = 2 w dw, w in [0, sqrt(xi)]
    static const double gx[8] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290,
                                 -0.1834346424956498, 0.1834346424956498,  0.5255324099163290,
                                 0.7966664774136267,  0.9602898564975363};
    static const double gw[8] = {0.1012285362903763, 0.2223810344533745, 0.3137066458778873,
                                 0.3626837833783620, 0.3626837833783620, 0.3137066458778873,
                                 0.2223810344533745, 0.1012285362903763};
    const double W = std::sqrt(xi);
    const int panels = 400;
    const double h = W / panels;
    const double se = std::sqrt(eta);
    double sum = 0.0;
    for (int p = 0; p < panels; ++p) {
        const double mid = (p + 0.5) * h;
        for (int g = 0; g < 8; ++g) {
            const double w = mid + 0.5 * h * gx[g];
            const double y = 2.0 * se * w;
            sum += gw[g] * 0.5 * h * 2.0 * w * std::exp(-(se - w) * (se - w)) * bessel_i0_scaled(y);
        }
    }
    return 1.0 - sum;
}

}  // namespace bed
