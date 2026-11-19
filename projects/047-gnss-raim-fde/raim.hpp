// Astron Built this AP
// Project 047 | Navigation | GNSS snapshot positioning with RAIM fault detection and exclusion
// Author: Aditya (aditya-astron)
// Date: 2026-11-19 IST
//
// Single-epoch GNSS position fix and receiver autonomous integrity monitoring.
//
//   rho_i = | R3(we tau_i) r_s,i - r | + b + e_i,   e_i ~ N(0, sigma^2)
//
// solved by Gauss-Newton for x = (r, b). The residual test uses
//   SSE / sigma^2 ~ chi2(N-4)            fault free
//   SSE / sigma^2 ~ ncchi2(N-4, lambda)  bias B on satellite i,
//                                        lambda = B^2 S_ii / sigma^2,
// with S = I - H (H^T H)^-1 H^T. A bias that gives noncentrality lambda moves
// the horizontal solution by slope_i * sigma * sqrt(lambda), where
//   slope_i = sqrt(A_E,i^2 + A_N,i^2) / sqrt(S_ii),  A = (H^T H)^-1 H^T.
// HPL = max_i slope_i * sigma * sqrt(lambda_det), lambda_det being the
// noncentrality at which the test misses with probability Pmd.
#pragma once

#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace raim {

constexpr double kPi = 3.14159265358979323846;
constexpr double kC = 299792458.0;            // m/s
constexpr double kOmegaE = 7.2921151467e-5;   // rad/s, WGS84 / IS-GPS-200
constexpr double kMuE = 3.986005e14;          // m^3/s^2, IS-GPS-200
constexpr double kWgsA = 6378137.0;
constexpr double kWgsF = 1.0 / 298.257223563;

using Vec3 = std::array<double, 3>;
using Mat = std::vector<std::vector<double>>;

inline double deg(double r) { return r * 180.0 / kPi; }
inline double rad(double d) { return d * kPi / 180.0; }
inline double norm(const Vec3& v) { return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]); }
inline Vec3 sub(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }

// ---------------------------------------------------------------- frames

Vec3 geodetic_to_ecef(double lat, double lon, double h);
// Rows are the E, N, U unit vectors expressed in ECEF.
std::array<Vec3, 3> enu_basis(double lat, double lon);
// Rotation about z by angle a (frame rotation, as used for the Sagnac term).
Vec3 rot_z(const Vec3& v, double a);

// ------------------------------------------------------------ linear algebra

Mat matmul(const Mat& a, const Mat& b);
Mat transpose(const Mat& a);
Mat inverse(Mat a);  // Gauss-Jordan with partial pivoting; throws if singular

// ------------------------------------------------------------- statistics

double gamma_p(double a, double x);               // regularised lower incomplete gamma
double chi2_cdf(double k, double x);
double chi2_inv(double k, double p);              // quantile
double ncchi2_cdf(double k, double lambda, double x);
// Noncentrality at which P(ncchi2(k, lambda) < t) = pmd.
double ncchi2_lambda_for_pmd(double k, double t, double pmd);

// ------------------------------------------------------------- positioning

struct Fix {
    Vec3 r{};            // ECEF position, m
    double b = 0.0;      // receiver clock bias, m
    int iterations = 0;
    std::vector<double> residuals;
    Mat H;               // ECEF geometry matrix at the solution, rows [-u^T 1]
};

// Gauss-Newton snapshot fix. sats are ECEF positions at transmit time,
// rho pseudoranges. With sagnac=true each satellite is rotated by we*tau
// into the ECEF frame at receive time.
Fix solve(const std::vector<Vec3>& sats, const std::vector<double>& rho, bool sagnac = true,
          Vec3 r0 = {0, 0, 0}, double b0 = 0.0, double tol = 1e-6, int max_iter = 20);

// Error-free pseudorange for a receiver at r with clock bias b.
double true_pseudorange(const Vec3& sat, const Vec3& r, double b, bool sagnac = true);

// Geometry matrix in ENU + clock from azimuth/elevation (radians).
Mat geometry_enu(const std::vector<double>& az, const std::vector<double>& el);
// Convert an ECEF geometry matrix to ENU columns at (lat, lon).
Mat to_enu(const Mat& H_ecef, double lat, double lon);

struct Dop {
    double gdop, pdop, hdop, vdop, tdop;
};
Dop dop(const Mat& H_enu);

// --------------------------------------------------------------- integrity

struct RaimConfig {
    double sigma = 5.0;    // pseudorange 1-sigma, m
    double pfa = 1e-5;     // per test
    double pmd = 1e-3;
};

struct RaimGeometry {
    int n = 0, dof = 0;
    double threshold = 0.0;     // on SSE / sigma^2
    double lambda_det = 0.0;
    std::vector<double> hslope, vslope, s_diag;
    std::vector<double> a_east, a_north, a_up;
    double hpl = INFINITY, vpl = INFINITY;
    int worst_h = -1;
    bool available = false;     // dof >= 1
};

RaimGeometry raim_geometry(const Mat& H_enu, const RaimConfig& cfg);

struct FdeResult {
    bool alarm = false;      // fault detected on the full set
    int excluded = -1;       // index removed, -1 if none
    bool recovered = false;  // a subset passed after exclusion
    double test = 0.0;       // SSE / sigma^2 on the full set
};

// Residual test, then exclusion of the satellite whose removal gives the
// smallest passing subset statistic (needs N >= 6 to be meaningful).
FdeResult detect_and_exclude(const std::vector<Vec3>& sats, const std::vector<double>& rho,
                             const RaimConfig& cfg, bool sagnac = true);

// --------------------------------------------------------------- constellation

struct Slot {
    double raan;   // rad, at t = 0
    double u0;     // argument of latitude at t = 0, rad
};

// Six planes at 55 deg, four unevenly spaced slots each, GPS-like layout.
std::vector<Slot> gps_like_24();
// ECEF position of a circular-orbit satellite at time t (GMST(0) = 0).
Vec3 sat_ecef(const Slot& s, double t, double incl = rad(55.0), double a = 26559.7e3);
void az_el(const Vec3& r_rx, double lat, double lon, const Vec3& r_sat, double& az, double& el);

}  // namespace raim
