// Astron Built this AP
// Project 015 | Electronics and circuit analysis | Damped LC input filter vs. constant-power-load instability
// Author: Aditya (aditya-astron)
// Date: 2026-11-05 IST
#include "filter.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {
constexpr double kPi = 3.14159265358979323846;
const double kGolden = (std::sqrt(5.0) - 1.0) / 2.0;

// Dense complex Gaussian elimination with partial pivoting. The matrices here
// are a handful of nodes, so nothing smarter is warranted.
std::vector<cplx> solve(std::vector<std::vector<cplx>> A, std::vector<cplx> b) {
    const size_t n = b.size();
    for (size_t k = 0; k < n; ++k) {
        size_t p = k;
        for (size_t r = k + 1; r < n; ++r)
            if (std::abs(A[r][k]) > std::abs(A[p][k])) p = r;
        if (std::abs(A[p][k]) == 0.0) throw std::runtime_error("singular nodal matrix");
        std::swap(A[k], A[p]);
        std::swap(b[k], b[p]);
        for (size_t r = k + 1; r < n; ++r) {
            const cplx m = A[r][k] / A[k][k];
            for (size_t c = k; c < n; ++c) A[r][c] -= m * A[k][c];
            b[r] -= m * b[k];
        }
    }
    std::vector<cplx> x(n);
    for (size_t k = n; k-- > 0;) {
        cplx s = b[k];
        for (size_t c = k + 1; c < n; ++c) s -= A[k][c] * x[c];
        x[k] = s / A[k][k];
    }
    return x;
}

template <class F>
double golden_max(F f, double a, double b, int iters = 100) {
    double c = b - kGolden * (b - a), d = a + kGolden * (b - a);
    double fc = f(c), fd = f(d);
    for (int k = 0; k < iters; ++k) {
        if (fc > fd) { b = d; d = c; fd = fc; c = b - kGolden * (b - a); fc = f(c); }
        else         { a = c; c = d; fc = fd; d = a + kGolden * (b - a); fd = f(d); }
    }
    return 0.5 * (a + b);
}
}  // namespace

cplx Netlist::impedance(int node, double f_hz) const {
    const cplx jw(0.0, 2.0 * kPi * f_hz);
    std::vector<std::vector<cplx>> Y(n_, std::vector<cplx>(n_, 0.0));
    for (const auto& e : els_) {
        cplx y;
        switch (e.kind) {
            case 'R': y = 1.0 / e.value; break;
            case 'L': y = 1.0 / (jw * e.value); break;
            default:  y = jw * e.value; break;
        }
        // Standard nodal stamp; ground rows/columns are dropped.
        const int a = e.a - 1, b = e.b - 1;
        if (a >= 0) Y[a][a] += y;
        if (b >= 0) Y[b][b] += y;
        if (a >= 0 && b >= 0) { Y[a][b] -= y; Y[b][a] -= y; }
    }
    std::vector<cplx> I(n_, 0.0);
    I[node - 1] = 1.0;
    return solve(Y, I)[node - 1];
}

Netlist build_netlist(const FilterDesign& d) {
    // 1: filter output / converter input. The shorted source is node 0. A zero
    // winding resistance cannot be stamped as a conductance, so the L then
    // goes straight to ground instead of through an internal node.
    int next = 2;
    const int nl = d.RL > 0.0 ? next++ : 0;
    const int nd = d.Cd > 0.0 ? next++ : -1;
    Netlist net(next - 1);
    if (d.RL > 0.0) net.R(0, nl, d.RL);
    net.L(nl, 1, d.L);
    net.C(1, 0, d.C);
    if (d.Cd > 0.0) {
        net.R(1, nd, d.Rd);
        net.C(nd, 0, d.Cd);
    }
    return net;
}

cplx zout_closed_form(const FilterDesign& d, double f_hz) {
    const cplx s(0.0, 2.0 * kPi * f_hz);
    cplx y = 1.0 / (d.RL + s * d.L) + s * d.C;
    if (d.Cd > 0.0) y += 1.0 / (d.Rd + 1.0 / (s * d.Cd));
    return 1.0 / y;
}

double peak_zout(const FilterDesign& d, double f_lo, double f_hi, double* f_peak) {
    const Netlist net = build_netlist(d);
    auto mag = [&](double lf) { return std::abs(net.impedance(1, std::pow(10.0, lf))); };
    const int N = 2000;
    const double a = std::log10(f_lo), b = std::log10(f_hi), h = (b - a) / N;
    int kmax = 0;
    double zmax = 0.0;
    for (int k = 0; k <= N; ++k) {
        const double z = mag(a + k * h);
        if (z > zmax) { zmax = z; kmax = k; }
    }
    const double lf = golden_max(mag, a + std::max(kmax - 1, 0) * h, a + std::min(kmax + 1, N) * h);
    if (f_peak) *f_peak = std::pow(10.0, lf);
    return mag(lf);
}

double r0(const FilterDesign& d) { return std::sqrt(d.L / d.C); }

double rd_optimal(double R0, double n) {
    return R0 * std::sqrt((2.0 + n) * (4.0 + 3.0 * n) / (2.0 * n * n * (4.0 + n)));
}

// At w = w0 sqrt(2/(2+n)) |Zout| is the same for Rd = 0 and Rd = inf, hence
// for every Rd; the optimum puts the peak exactly on this fixed point.
double zpeak_optimal(double R0, double n) { return R0 * std::sqrt(2.0 * (2.0 + n)) / n; }

double fpeak_optimal(double f0, double n) { return f0 * std::sqrt(2.0 / (2.0 + n)); }

double rd_minimax(FilterDesign d, double f_lo, double f_hi) {
    // Too little Rd leaves the Cd resonance undamped, too much leaves the C
    // resonance undamped; peak |Zout| vs Rd is unimodal with a minimum between.
    auto neg_peak = [&](double lr) {
        d.Rd = std::pow(10.0, lr);
        return -peak_zout(d, f_lo, f_hi);
    };
    const double R0 = r0(d);
    return std::pow(10.0, golden_max(neg_peak, std::log10(0.05 * R0), std::log10(5.0 * R0), 60));
}

CplOperatingPoint operating_point(const FilterDesign& d, double vin, double p) {
    // DC: Vin - RL*P/v = v  ->  v^2 - Vin v + RL P = 0, take the high-voltage root.
    const double disc = vin * vin - 4.0 * d.RL * p;
    if (disc < 0.0) return {false, 0.0, 0.0, 0.0};
    const double v = 0.5 * (vin + std::sqrt(disc));
    return {true, v, p / v, v * v / p};
}

bool is_stable(const FilterDesign& d, double vin, double p) {
    const auto op = operating_point(d, vin, p);
    if (!op.exists) return false;
    // States x = [iL, vC, vCd]. Linearised CPL current: di = -(1/Rn) dv.
    //   L  diL/dt  = -RL iL - vC
    //   C  dvC/dt  = iL + vC/Rn - (vC - vCd)/Rd
    //   Cd dvCd/dt = (vC - vCd)/Rd
    const double gn = 1.0 / op.rn;
    if (d.Cd <= 0.0) {
        // s^2 + a1 s + a0
        const double a1 = d.RL / d.L - gn / d.C;
        const double a0 = (1.0 - d.RL * gn) / (d.L * d.C);
        return a1 > 0.0 && a0 > 0.0;
    }
    const double gd = 1.0 / d.Rd;
    const double A[3][3] = {
        {-d.RL / d.L, -1.0 / d.L, 0.0},
        {1.0 / d.C, (gn - gd) / d.C, gd / d.C},
        {0.0, gd / d.Cd, -gd / d.Cd},
    };
    // det(sI - A) = s^3 + a2 s^2 + a1 s + a0
    const double tr = A[0][0] + A[1][1] + A[2][2];
    const double m2 = A[0][0] * A[1][1] - A[0][1] * A[1][0] +
                      A[0][0] * A[2][2] - A[0][2] * A[2][0] +
                      A[1][1] * A[2][2] - A[1][2] * A[2][1];
    const double det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
                       A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
                       A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
    const double a2 = -tr, a1 = m2, a0 = -det;
    return a2 > 0.0 && a1 > 0.0 && a0 > 0.0 && a2 * a1 > a0;
}

double max_stable_power(const FilterDesign& d, double vin, double p_hi) {
    double lo = 1e-6, hi = p_hi;
    if (is_stable(d, vin, hi)) return hi;
    for (int k = 0; k < 200 && hi - lo > 1e-9 * hi; ++k) {
        const double mid = 0.5 * (lo + hi);
        (is_stable(d, vin, mid) ? lo : hi) = mid;
    }
    return lo;
}

double simulate_growth(const FilterDesign& d, double vin, double p, double dv,
                       double t_end, double dt) {
    const auto op = operating_point(d, vin, p);
    if (!op.exists) return std::numeric_limits<double>::infinity();
    const bool damped = d.Cd > 0.0;
    auto f = [&](const double x[3], double out[3]) {
        const double id = damped ? (x[1] - x[2]) / d.Rd : 0.0;
        out[0] = (vin - d.RL * x[0] - x[1]) / d.L;
        out[1] = (x[0] - p / x[1] - id) / d.C;
        out[2] = damped ? id / d.Cd : 0.0;
    };
    double x[3] = {op.i, op.v + dv, op.v};
    const long steps = static_cast<long>(t_end / dt);
    const long tail = steps - steps / 10;
    double peak = 0.0;
    for (long k = 0; k < steps; ++k) {
        double k1[3], k2[3], k3[3], k4[3], t[3];
        f(x, k1);
        for (int j = 0; j < 3; ++j) t[j] = x[j] + 0.5 * dt * k1[j];
        f(t, k2);
        for (int j = 0; j < 3; ++j) t[j] = x[j] + 0.5 * dt * k2[j];
        f(t, k3);
        for (int j = 0; j < 3; ++j) t[j] = x[j] + dt * k3[j];
        f(t, k4);
        for (int j = 0; j < 3; ++j) x[j] += dt / 6.0 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
        // A real converter drops out under UVLO long before the CPL model
        // goes singular at v = 0.
        if (x[1] < 0.5 * vin) return std::numeric_limits<double>::infinity();
        if (k >= tail) peak = std::max(peak, std::abs(x[1] - op.v));
    }
    return peak / dv;
}
