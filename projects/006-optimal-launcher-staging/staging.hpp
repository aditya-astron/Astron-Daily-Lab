// Astron Built this AP
// Project 006 | Rocket propulsion and staging | Optimal mass staging of a multistage launcher
// Author: Aditya (aditya-astron)
// Date: 2026-11-02 IST
//
// Restricted staging model for serial (tandem) stages:
//   structural ratio   eps_i = m_inert_i / (m_inert_i + m_prop_i)
//   mass ratio         n_i   = m_ignition_i / m_burnout_i
//   delta-v            dv    = sum c_i ln n_i,  c_i = g0 * Isp_i
// Everything above stage i (upper stages + payload) is that stage's payload.

#pragma once

#include <string>
#include <vector>

namespace staging {

constexpr double kG0 = 9.80665;  // m/s^2, standard gravity (defines Isp in s)

struct Stage {
    std::string name;
    double isp_s;  // effective specific impulse over the burn
    double eps;    // structural ratio, 0 < eps < 1
    double c() const { return kG0 * isp_s; }
};

struct StageMass {
    double ignition;   // stage + everything above it
    double inert;
    double propellant;
    double dv;         // velocity gained in this stage's burn
    double n;          // mass ratio
};

struct Solution {
    bool ok = false;
    std::string error;
    double eta = 0.0;          // Lagrange multiplier [s/m]
    double gross = 0.0;        // lift-off mass [kg]
    double payload = 0.0;
    std::vector<StageMass> stages;  // index 0 = first (bottom) stage
};

// Upper bound on dv for a stack: every stage burning to n_i -> 1/eps_i
// (zero payload limit).
double max_delta_v(const std::vector<Stage>& s);

// Build stage masses for a prescribed set of mass ratios, top-down from the
// payload. Fails if any n_i is outside (1, 1/eps_i).
Solution evaluate(const std::vector<Stage>& s, const std::vector<double>& n,
                  double payload);

// Minimum lift-off mass for a required dv (Lagrange multiplier method).
Solution optimize(const std::vector<Stage>& s, double dv, double payload);

// Equal dv per stage, for comparison against the optimum.
Solution equal_split(const std::vector<Stage>& s, double dv, double payload);

}  // namespace staging
