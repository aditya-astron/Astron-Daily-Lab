// Astron Built this AP
// Project 015 | Electronics and circuit analysis | Damped LC input filter vs. constant-power-load instability
// Author: Aditya (aditya-astron)
// Date: 2026-11-05 IST
#pragma once
#include <complex>
#include <vector>

using cplx = std::complex<double>;

// Small-signal AC nodal solver. Node 0 is ground. Only passive two-terminal
// elements are needed: the input source is an ideal voltage source, so for
// output impedance it is replaced by a short and the "Vin" node becomes ground.
class Netlist {
public:
    explicit Netlist(int nodes) : n_(nodes) {}
    void R(int a, int b, double ohm) { els_.push_back({'R', a, b, ohm}); }
    void L(int a, int b, double henry) { els_.push_back({'L', a, b, henry}); }
    void C(int a, int b, double farad) { els_.push_back({'C', a, b, farad}); }

    // Driving-point impedance at `node`: inject 1 A, solve Y v = i, read v.
    cplx impedance(int node, double f_hz) const;

private:
    struct Element { char kind; int a, b; double value; };
    int n_;
    std::vector<Element> els_;
};

// Single-stage LC filter with an Rd-Cd parallel damping leg (Middlebrook /
// Erickson topology). RL is the inductor winding resistance. Cd = 0 means
// no damping leg.
struct FilterDesign {
    double L, C, RL, Rd, Cd;
};

Netlist build_netlist(const FilterDesign& d);         // output node = 1
cplx zout_closed_form(const FilterDesign& d, double f_hz);

// Peak |Zout| over [f_lo, f_hi], log sweep then golden-section refinement.
double peak_zout(const FilterDesign& d, double f_lo, double f_hi, double* f_peak = nullptr);

// Erickson, Fundamentals of Power Electronics, 2nd ed., sec. 10.4.1:
// with n = Cd/C and R0 = sqrt(L/C), the Rd that minimises the peak output
// impedance and the resulting peak, assuming a lossless inductor.
double r0(const FilterDesign& d);
double rd_optimal(double R0, double n);
double zpeak_optimal(double R0, double n);
double fpeak_optimal(double f0, double n);

// Rd that minimises the peak |Zout| found numerically (golden section).
double rd_minimax(FilterDesign d, double f_lo, double f_hi);

// Converter modelled as an ideal constant-power load P fed from Vin through
// the filter. Small-signal it is a negative resistance -Rn = -V^2/P.
struct CplOperatingPoint {
    bool exists;   // false if P exceeds the max-power-transfer limit Vin^2/(4 RL)
    double v, i, rn;
};
CplOperatingPoint operating_point(const FilterDesign& d, double vin, double p);

// Routh-Hurwitz on the linearised filter + CPL state equations.
bool is_stable(const FilterDesign& d, double vin, double p);

// Highest load power that is still small-signal stable (bisection on is_stable).
double max_stable_power(const FilterDesign& d, double vin, double p_hi);

// Nonlinear RK4 time simulation: start at the operating point, kick the
// capacitor voltage by dv, run for t_end. Returns the ratio of the peak
// |v - v_eq| over the last 10% of the run to the initial kick. A collapse of
// the bus below vin/2 returns +inf.
double simulate_growth(const FilterDesign& d, double vin, double p, double dv,
                       double t_end, double dt);
