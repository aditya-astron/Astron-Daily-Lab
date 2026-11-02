// Astron Built this AP
// Project 006 | Rocket propulsion and staging | Optimal mass staging of a multistage launcher
// Author: Aditya (aditya-astron)
// Date: 2026-11-02 IST
//
// Optimality condition. Per stage, with lambda_i = payload_i / (stage_i):
//   m_ign_i / m_pl_i = n_i (1 - eps_i) / (1 - eps_i n_i)
// so  ln(m0 / m_PL) = sum [ ln n_i + ln(1 - eps_i) - ln(1 - eps_i n_i) ].
// Minimising this subject to sum c_i ln n_i = dv, with x_i = ln n_i, gives
//   1 / (1 - eps_i n_i) = eta c_i   =>   n_i = (eta c_i - 1) / (eta c_i eps_i)
// and eta is the root of  f(eta) = sum c_i ln n_i(eta) - dv.
// f is strictly increasing in eta, so a bracketed bisection is robust.
// (Curtis, Orbital Mechanics for Engineering Students, sec. 11.6.)

#include "staging.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace staging {

double max_delta_v(const std::vector<Stage>& s) {
    double v = 0.0;
    for (const auto& st : s) v += st.c() * std::log(1.0 / st.eps);
    return v;
}

Solution evaluate(const std::vector<Stage>& s, const std::vector<double>& n,
                  double payload) {
    Solution sol;
    sol.payload = payload;
    if (n.size() != s.size()) {
        sol.error = "mass ratio count does not match stage count";
        return sol;
    }
    sol.stages.resize(s.size());
    double above = payload;
    for (std::size_t k = s.size(); k-- > 0;) {
        const double e = s[k].eps;
        if (!(n[k] > 1.0) || !(n[k] * e < 1.0)) {
            sol.error = "stage " + s[k].name + ": mass ratio outside (1, 1/eps)";
            return sol;
        }
        const double ign = above * n[k] * (1.0 - e) / (1.0 - e * n[k]);
        const double stage = ign - above;
        sol.stages[k] = {ign, e * stage, (1.0 - e) * stage,
                         s[k].c() * std::log(n[k]), n[k]};
        above = ign;
    }
    sol.gross = above;
    sol.ok = true;
    return sol;
}

Solution optimize(const std::vector<Stage>& s, double dv, double payload) {
    Solution bad;
    bad.payload = payload;
    if (s.empty() || dv <= 0.0 || payload <= 0.0) {
        bad.error = "need at least one stage, dv > 0 and payload > 0";
        return bad;
    }
    if (dv >= max_delta_v(s)) {
        bad.error = "required dv exceeds the zero-payload limit of the stack";
        return bad;
    }

    auto ratios = [&](double eta) {
        std::vector<double> n(s.size());
        for (std::size_t i = 0; i < s.size(); ++i) {
            const double ec = eta * s[i].c();
            n[i] = (ec - 1.0) / (ec * s[i].eps);
        }
        return n;
    };
    auto f = [&](double eta) {
        double v = 0.0;
        for (std::size_t i = 0; i < s.size(); ++i) {
            const double ec = eta * s[i].c();
            v += s[i].c() * std::log((ec - 1.0) / (ec * s[i].eps));
        }
        return v - dv;
    };

    // n_i > 1 requires eta > 1 / (c_i (1 - eps_i)) for every stage.
    double lo = 0.0;
    for (const auto& st : s) lo = std::max(lo, 1.0 / (st.c() * (1.0 - st.eps)));
    lo *= 1.0 + 1e-12;

    // At the lower bound the weakest stage contributes nothing. If the others
    // already exceed dv there, the optimum wants that stage removed.
    if (f(lo) > 0.0) {
        bad.error = "a stage would have n <= 1 at the optimum; drop a stage";
        return bad;
    }
    double hi = 2.0 * lo;
    while (f(hi) < 0.0) hi *= 2.0;  // terminates: f(inf) = max_dv - dv > 0

    for (int it = 0; it < 200 && (hi - lo) > 1e-15 * hi; ++it) {
        const double mid = 0.5 * (lo + hi);
        (f(mid) < 0.0 ? lo : hi) = mid;
    }
    const double eta = 0.5 * (lo + hi);

    Solution sol = evaluate(s, ratios(eta), payload);
    sol.eta = eta;
    return sol;
}

Solution equal_split(const std::vector<Stage>& s, double dv, double payload) {
    std::vector<double> n(s.size());
    const double per = dv / static_cast<double>(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) n[i] = std::exp(per / s[i].c());
    return evaluate(s, n, payload);
}

}  // namespace staging
