// Astron Built this AP
// Project 018 | Numerical methods | Symplectic vs. RK4 integration of long-arc Keplerian orbits
// Author: Aditya (aditya-astron)
// Date: 2026-11-06 IST
#include "integrators.hpp"

#include <stdexcept>

namespace orb {

State step_rk4(const State& s, double h, double mu) {
    // y' = (v, a(r))
    const Vec3 k1r = s.v;
    const Vec3 k1v = gravity(s.r, mu);
    const Vec3 k2r = s.v + (0.5 * h) * k1v;
    const Vec3 k2v = gravity(s.r + (0.5 * h) * k1r, mu);
    const Vec3 k3r = s.v + (0.5 * h) * k2v;
    const Vec3 k3v = gravity(s.r + (0.5 * h) * k2r, mu);
    const Vec3 k4r = s.v + h * k3v;
    const Vec3 k4v = gravity(s.r + h * k3r, mu);
    const double c = h / 6.0;
    return {s.r + c * (k1r + 2.0 * k2r + 2.0 * k3r + k4r),
            s.v + c * (k1v + 2.0 * k2v + 2.0 * k3v + k4v)};
}

State step_leapfrog(const State& s, double h, double mu) {
    // Each sub-map is an exact flow of part of H = T(v) + V(r), so the
    // composition is symplectic. Kicks are radial and drifts move r along v,
    // so r x v is preserved exactly (to round-off) for any central force.
    const Vec3 v_half = s.v + (0.5 * h) * gravity(s.r, mu);
    const Vec3 r_new = s.r + h * v_half;
    const Vec3 v_new = v_half + (0.5 * h) * gravity(r_new, mu);
    return {r_new, v_new};
}

State step_yoshida4(const State& s, double h, double mu) {
    static const double cbrt2 = std::cbrt(2.0);
    static const double w1 = 1.0 / (2.0 - cbrt2);
    static const double w0 = -cbrt2 / (2.0 - cbrt2);
    State x = step_leapfrog(s, w1 * h, mu);
    x = step_leapfrog(x, w0 * h, mu);
    return step_leapfrog(x, w1 * h, mu);
}

State kepler_propagate(const State& s0, double dt, double mu) {
    const double r0 = norm(s0.r);
    const double v0sq = dot(s0.v, s0.v);
    const double a = 1.0 / (2.0 / r0 - v0sq / mu);  // vis-viva
    if (!(a > 0.0)) throw std::domain_error("kepler_propagate: orbit is not elliptic");

    const double sqa = std::sqrt(a);
    const double n = std::sqrt(mu / (a * a * a));
    const double sigma0 = dot(s0.r, s0.v) / std::sqrt(mu);
    const double ecosE0 = 1.0 - r0 / a;      // e cos E0
    const double esinE0 = sigma0 / sqa;      // e sin E0

    // Reduce the mean-anomaly change to one revolution, solve, then add the
    // whole revolutions back. Keeps Newton well conditioned on long arcs.
    const double two_pi = 2.0 * M_PI;
    const double dM_full = n * dt;
    const double revs = std::floor(dM_full / two_pi);
    const double dM = dM_full - revs * two_pi;

    double dE = dM;  // starting guess, fine for e < ~0.9 with the safeguard below
    for (int it = 0; it < 60; ++it) {
        const double sd = std::sin(dE), cd = std::cos(dE);
        const double F = dE - ecosE0 * sd + esinE0 * (1.0 - cd) - dM;
        const double dF = 1.0 - ecosE0 * cd + esinE0 * sd;  // = r/a > 0
        double step = F / dF;
        if (std::fabs(step) > 1.0) step = step > 0 ? 1.0 : -1.0;
        dE -= step;
        if (std::fabs(step) < 1e-15) break;
    }

    const double sd = std::sin(dE), cd = std::cos(dE);
    const double r = a + (r0 - a) * cd + sigma0 * sqa * sd;
    const double dE_total = dE + revs * two_pi;

    const double f = 1.0 - a / r0 * (1.0 - cd);
    const double g = dt - std::sqrt(a * a * a / mu) * (dE_total - sd);
    const double fdot = -std::sqrt(mu * a) / (r * r0) * sd;
    const double gdot = 1.0 - a / r * (1.0 - cd);

    return {f * s0.r + g * s0.v, fdot * s0.r + gdot * s0.v};
}

State from_elements(double a, double e, double inc, double raan, double argp,
                    double mean_anom, double mu) {
    if (!(a > 0.0) || !(e >= 0.0 && e < 1.0)) throw std::domain_error("from_elements: elliptic only");
    double E = mean_anom;
    for (int it = 0; it < 60; ++it) {
        const double dE = (E - e * std::sin(E) - mean_anom) / (1.0 - e * std::cos(E));
        E -= dE;
        if (std::fabs(dE) < 1e-15) break;
    }
    const double b = a * std::sqrt(1.0 - e * e);
    const double r = a * (1.0 - e * std::cos(E));
    const double Edot = std::sqrt(mu / (a * a * a)) * a / r;
    // Perifocal (P toward periapsis, Q 90 deg ahead in the orbit plane).
    const double xp = a * (std::cos(E) - e), yp = b * std::sin(E);
    const double vxp = -a * std::sin(E) * Edot, vyp = b * std::cos(E) * Edot;

    const double cO = std::cos(raan), sO = std::sin(raan);
    const double ci = std::cos(inc), si = std::sin(inc);
    const double cw = std::cos(argp), sw = std::sin(argp);
    const Vec3 P = {cO * cw - sO * sw * ci, sO * cw + cO * sw * ci, sw * si};
    const Vec3 Q = {-cO * sw - sO * cw * ci, -sO * sw + cO * cw * ci, cw * si};
    return {xp * P + yp * Q, vxp * P + vyp * Q};
}

const char* method_name(Method m) {
    switch (m) {
        case Method::RK4: return "RK4";
        case Method::Leapfrog: return "Leapfrog";
        case Method::Yoshida4: return "Yoshida4";
    }
    return "?";
}

int force_evals_per_step(Method m) {
    // Leapfrog's closing kick force is reused as the next opening kick when
    // stepping repeatedly, so its steady-state cost is one evaluation.
    switch (m) {
        case Method::RK4: return 4;
        case Method::Leapfrog: return 1;
        case Method::Yoshida4: return 3;
    }
    return 0;
}

State step(Method m, const State& s, double h, double mu) {
    switch (m) {
        case Method::RK4: return step_rk4(s, h, mu);
        case Method::Leapfrog: return step_leapfrog(s, h, mu);
        case Method::Yoshida4: return step_yoshida4(s, h, mu);
    }
    throw std::invalid_argument("unknown method");
}

}  // namespace orb
