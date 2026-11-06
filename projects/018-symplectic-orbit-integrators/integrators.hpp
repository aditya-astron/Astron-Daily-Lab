// Astron Built this AP
// Project 018 | Numerical methods | Symplectic vs. RK4 integration of long-arc Keplerian orbits
// Author: Aditya (aditya-astron)
// Date: 2026-11-06 IST
#pragma once

#include <array>
#include <cmath>

namespace orb {

using Vec3 = std::array<double, 3>;

inline Vec3 operator+(const Vec3& a, const Vec3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
inline Vec3 operator-(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
inline Vec3 operator*(double s, const Vec3& a) { return {s * a[0], s * a[1], s * a[2]}; }
inline double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
inline double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

struct State {
    Vec3 r;  // position [km]
    Vec3 v;  // velocity [km/s]
};

// Two-body acceleration a = -mu r / |r|^3.
inline Vec3 gravity(const Vec3& r, double mu) {
    const double rn = norm(r);
    return (-mu / (rn * rn * rn)) * r;
}

// Specific orbital energy E = v^2/2 - mu/r and angular momentum h = r x v.
inline double energy(const State& s, double mu) { return 0.5 * dot(s.v, s.v) - mu / norm(s.r); }
inline Vec3 ang_mom(const State& s) { return cross(s.r, s.v); }

// Classical 4th-order Runge-Kutta. Four force evaluations per step.
// Not symplectic: the energy error accumulates secularly.
State step_rk4(const State& s, double h, double mu);

// Kick-drift-kick leapfrog (velocity Verlet). 2nd order, symplectic,
// time-reversible; one force evaluation per step when the last kick is reused.
State step_leapfrog(const State& s, double h, double mu);

// Yoshida (1990) 4th-order composition: three leapfrog substeps with
//   w1 = 1 / (2 - 2^(1/3)),  w0 = -2^(1/3) / (2 - 2^(1/3)),  w1 + w0 + w1 = 1.
// The middle substep runs backward in time; the symmetric composition
// cancels the h^3 error term of leapfrog. Three force evaluations per step.
State step_yoshida4(const State& s, double h, double mu);

// Exact two-body propagation using Lagrange f and g coefficients in terms of
// the eccentric anomaly change dE (elliptic orbits only).
// Solves  n dt = dE - (1 - r0/a) sin dE + sigma0/sqrt(a) (1 - cos dE),
// sigma0 = r0.v0 / sqrt(mu), by Newton iteration.
State kepler_propagate(const State& s0, double dt, double mu);

// Cartesian state from classical elements (angles in radians), via the
// perifocal frame and the 3-1-3 rotation R3(-raan) R1(-i) R3(-argp).
State from_elements(double a, double e, double inc, double raan, double argp,
                    double mean_anom, double mu);

enum class Method { RK4, Leapfrog, Yoshida4 };

const char* method_name(Method m);
int force_evals_per_step(Method m);
State step(Method m, const State& s, double h, double mu);

}  // namespace orb
