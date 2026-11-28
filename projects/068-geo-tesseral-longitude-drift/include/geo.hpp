// Astron Built this AP
// Project 068 | Numerical methods / Orbital mechanics | Spherical-harmonic gravity and GEO longitude drift
// Author: Aditya (aditya-astron)
// Date: 2026-11-28 IST
//
// Longitude dynamics of a geostationary satellite under the tesseral field.
//
// A tangential acceleration a_T changes the semi-major axis at da/dt = 2 a_T / n.
// The drift rate relative to the Earth is lam_dot = n - wE, and dn/da = -3n/(2a), so
//
//   lam_ddot = -3 a_T / a
//
// Equilibria sit where a_T(lam) = 0; they are stable where d(lam_ddot)/d(lam) < 0.
// The 3D propagator below integrates the full orbit in inertial space with the
// Earth rotating uniformly, and is used to check the one-dimensional reduction.

#pragma once
#include <functional>
#include <vector>

#include "gravity.hpp"

namespace geo {

using grav::Field;
using grav::Vec3;

constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = kPi / 180.0;
constexpr double kWE = 7.292115e-5;             // Earth rotation rate, rad/s
constexpr double kSidDay = 2.0 * kPi / kWE;     // 86164.09 s
constexpr double kDay = 86400.0;
constexpr double kYear = 365.25 * kDay;

inline double sync_radius(const Field& f) { return std::cbrt(f.GM / (kWE * kWE)); }

inline double wrap360(double deg) {
    deg = std::fmod(deg, 360.0);
    return deg < 0 ? deg + 360.0 : deg;
}

// Along-track (east) acceleration on the equator at longitude lam (rad). Taken
// straight from (1/r) dU/dlam: projecting the Cartesian vector would bury this
// 1e-8 m/s^2 signal under round-off from the 0.22 m/s^2 central term.
inline double along_track(const Field& f, double lam, double r) {
    return grav::accel_local(f, r, 0.0, lam)[2];
}

inline double lam_ddot(const Field& f, double lam, double r) { return -3.0 * along_track(f, lam, r) / r; }

struct Equilibrium { double lam_deg; bool stable; double curvature; };   // curvature = d(lam_ddot)/d(lam), 1/s^2

inline std::vector<Equilibrium> equilibria(const Field& f, double r) {
    std::vector<Equilibrium> out;
    const int K = 3600;
    auto g = [&](double d) { return along_track(f, d * kDeg, r); };
    double prev = g(0.0);
    for (int i = 1; i <= K; ++i) {
        double lo = (i - 1) * 360.0 / K, hi = i * 360.0 / K, cur = g(hi);
        if ((prev < 0) != (cur < 0)) {
            double flo = prev;
            for (int it = 0; it < 80; ++it) {
                double mid = 0.5 * (lo + hi), fm = g(mid);
                if ((fm < 0) == (flo < 0)) { lo = mid; flo = fm; } else { hi = mid; }
            }
            double l = 0.5 * (lo + hi), h = 1e-3;
            double k = (lam_ddot(f, (l + h) * kDeg, r) - lam_ddot(f, (l - h) * kDeg, r)) / (2 * h * kDeg);
            out.push_back({wrap360(l), k < 0, k});
        }
        prev = cur;
    }
    return out;
}

// Period of the 1D model lam_ddot = F(lam) released at rest from lam0,
// measured from successive upward crossings of the equilibrium.
inline double libration_period_1d(const Field& f, double r, double lam0, double lam_eq) {
    const double dt = 0.05 * kDay;
    double lam = lam0, w = 0.0, t = 0.0, last = -1.0, sum = 0.0;
    int n = 0;
    auto F = [&](double l) { return lam_ddot(f, l, r); };
    while (n < 4 && t < 1e5 * kDay) {
        double l0 = lam;
        double k1l = w, k1w = F(lam);
        double k2l = w + 0.5 * dt * k1w, k2w = F(lam + 0.5 * dt * k1l);
        double k3l = w + 0.5 * dt * k2w, k3w = F(lam + 0.5 * dt * k2l);
        double k4l = w + dt * k3w, k4w = F(lam + dt * k3l);
        lam += dt / 6 * (k1l + 2 * k2l + 2 * k3l + k4l);
        w += dt / 6 * (k1w + 2 * k2w + 2 * k3w + k4w);
        t += dt;
        if (l0 < lam_eq && lam >= lam_eq) {
            double tc = t - dt + dt * (lam_eq - l0) / (lam - l0);
            if (last >= 0) { sum += tc - last; ++n; }
            last = tc;
        }
    }
    return n ? sum / n : 0.0;
}

// ---------------------------------------------------------------- 3D model

struct State { Vec3 r, v; };

inline Vec3 accel_inertial(const Field& f, double t, const Vec3& x) {
    double th = kWE * t, c = std::cos(th), s = std::sin(th);
    Vec3 xe{c * x[0] + s * x[1], -s * x[0] + c * x[1], x[2]};
    Vec3 ae = grav::accel_spherical(f, xe);
    return {c * ae[0] - s * ae[1], s * ae[0] + c * ae[1], ae[2]};
}

inline void rk4_step(const Field& f, double t, double dt, State& y) {
    auto deriv = [&](double tt, const State& s) {
        return State{s.v, accel_inertial(f, tt, s.r)};
    };
    auto add = [](const State& a, const State& b, double h) {
        State o;
        for (int i = 0; i < 3; ++i) { o.r[i] = a.r[i] + h * b.r[i]; o.v[i] = a.v[i] + h * b.v[i]; }
        return o;
    };
    State k1 = deriv(t, y);
    State k2 = deriv(t + dt / 2, add(y, k1, dt / 2));
    State k3 = deriv(t + dt / 2, add(y, k2, dt / 2));
    State k4 = deriv(t + dt, add(y, k3, dt));
    for (int i = 0; i < 3; ++i) {
        y.r[i] += dt / 6 * (k1.r[i] + 2 * k2.r[i] + 2 * k3.r[i] + k4.r[i]);
        y.v[i] += dt / 6 * (k1.v[i] + 2 * k2.v[i] + 2 * k3.v[i] + k4.v[i]);
    }
}

// Earth-fixed longitude of an inertial position at time t.
inline double ecef_longitude(double t, const Vec3& x) { return std::atan2(x[1], x[0]) - kWE * t; }

// Circular equatorial start at lam0 with zero drift: pick r so the local radial
// gravity (all harmonics) supplies exactly wE^2 r, then v = wE r eastward.
inline State geostationary_start(const Field& f, double lam0) {
    auto resid = [&](double r) {
        Vec3 x{r * std::cos(lam0), r * std::sin(lam0), 0.0};
        Vec3 a = grav::accel_spherical(f, x);
        double ar = (a[0] * x[0] + a[1] * x[1]) / r;
        return -ar - kWE * kWE * r;
    };
    double lo = sync_radius(f) - 5e3, hi = sync_radius(f) + 5e3;
    for (int i = 0; i < 100; ++i) {
        double mid = 0.5 * (lo + hi);
        if ((resid(mid) > 0) == (resid(lo) > 0)) lo = mid; else hi = mid;
    }
    double r = 0.5 * (lo + hi), v = kWE * r;
    return {{r * std::cos(lam0), r * std::sin(lam0), 0.0}, {-v * std::sin(lam0), v * std::cos(lam0), 0.0}};
}

// Sidereal-day running mean of the unwrapped Earth-fixed longitude. Steps per
// sidereal day is an integer, so the mean removes the daily eccentricity wobble
// exactly for a fixed orbit.
class DailyMean {
public:
    explicit DailyMean(int n) : buf_(n, 0.0) {}
    void push(double lam) {
        if (count_ > 0) {   // unwrap against the previous sample
            while (lam - last_ > kPi) lam -= 2 * kPi;
            while (lam - last_ < -kPi) lam += 2 * kPi;
        }
        last_ = lam;
        sum_ += lam - buf_[head_];
        buf_[head_] = lam;
        head_ = (head_ + 1) % buf_.size();
        ++count_;
    }
    bool ready() const { return count_ >= buf_.size(); }
    double mean() const { return sum_ / buf_.size(); }
private:
    std::vector<double> buf_;
    size_t head_ = 0, count_ = 0;
    double sum_ = 0.0, last_ = 0.0;
};

struct Libration3D { double period; double min_deg, max_deg; };

// Release at rest from lam0 and time the daily-mean longitude oscillation.
inline Libration3D libration_period_3d(const Field& f, double lam0, double lam_eq, int crossings = 3) {
    const int per_day = 720;
    const double dt = kSidDay / per_day;
    State y = geostationary_start(f, lam0);
    DailyMean dm(per_day);
    double t = 0.0, prev = 0.0, last = -1.0, sum = 0.0, lo = 1e9, hi = -1e9;
    int n = 0;
    bool have_prev = false;
    while (n < crossings && t < 20000 * kDay) {
        rk4_step(f, t, dt, y);
        t += dt;
        dm.push(ecef_longitude(t, y.r));
        if (!dm.ready()) continue;
        double m = dm.mean();
        // the mean is centred half a sidereal day back
        double tc_now = t - 0.5 * kSidDay;
        double rel = std::remainder(m - lam_eq, 2 * kPi);
        lo = std::min(lo, rel); hi = std::max(hi, rel);
        if (have_prev && prev < 0 && rel >= 0) {
            double tc = tc_now - dt + dt * (-prev) / (rel - prev);
            if (last >= 0) { sum += tc - last; ++n; }
            last = tc;
        }
        prev = rel;
        have_prev = true;
    }
    return {n ? sum / n : 0.0, lo / kDeg, hi / kDeg};
}

struct StationKeeping {
    int burns = 0;
    double dv_total = 0.0;        // m/s, all burns
    double dv_rate = 0.0;         // m/s per year, measured over whole cycles
    double max_excursion = 0.0;   // deg from slot centre, daily-mean longitude
    double cycle_days = 0.0;
};

// East-west box keeping with paired tangential burns 12 h apart (the pair
// leaves eccentricity unchanged to first order). lam_c and half_box in radians.
//
// Each cycle is a parabola in longitude: leaving the downhill edge with speed
// v = sqrt(2 A W) against the drift acceleration A carries the satellite
// exactly to the uphill edge (W = full box width) and back, in 2 v / A.
// Per cycle dv = 2 v a / 3, so the cost rate is a A / 3 = |a_T|, independent
// of the box width. The drift is
// estimated from the daily-mean longitude, so the burn sets the post-burn
// drift to -v instead of applying a fixed open-loop reversal; a fixed
// reversal lets the timing error of each cycle feed the next.
inline StationKeeping east_west_keeping(const Field& f, double lam_c, double half_box, double years) {
    const int per_day = 720;
    const double dt = kSidDay / per_day;
    const double r0 = sync_radius(f);
    const double A = lam_ddot(f, lam_c, r0);          // rad/s^2, sign = downhill direction
    const double sgn = A > 0 ? 1.0 : -1.0;
    const double v_target = std::sqrt(2.0 * std::fabs(A) * 2.0 * half_box);

    State y = geostationary_start(f, lam_c - sgn * half_box);
    DailyMean dm(per_day);
    std::vector<double> hist(per_day, 0.0);           // daily means, one day back
    StationKeeping out;
    double t = 0.0, first_burn = -1.0, last_burn = -1.0, dv_cycles = 0.0, half = 0.0;
    int pending_half = -1, lockout = 0, pairs = 0;
    size_t k = 0;
    auto burn = [&](double dv) {   // along +velocity when dv > 0
        double vn = std::sqrt(y.v[0] * y.v[0] + y.v[1] * y.v[1] + y.v[2] * y.v[2]);
        for (int i = 0; i < 3; ++i) y.v[i] += dv * y.v[i] / vn;
        out.dv_total += std::fabs(dv);
        ++out.burns;
    };
    while (t < years * kYear) {
        rk4_step(f, t, dt, y);
        t += dt;
        dm.push(ecef_longitude(t, y.r));
        if (pending_half > 0 && --pending_half == 0) burn(half);
        if (lockout > 0) --lockout;
        if (!dm.ready()) continue;

        double m = dm.mean();
        double old = hist[k % per_day];
        hist[k % per_day] = m;
        ++k;
        if (k <= size_t(per_day)) continue;            // need a day of means for the slope

        double rel = std::remainder(m - lam_c, 2 * kPi);
        out.max_excursion = std::max(out.max_excursion, std::fabs(rel) / kDeg);
        if (lockout > 0) continue;

        // Slope of the mean over the last day is the drift one day ago (the
        // mean is centred half a day back); push it to now with A.
        double rate = (m - old) / kSidDay + A * kSidDay;
        // Where the satellite will be at the centroid of the burn pair:
        // half a day of mean lag plus a quarter day to the pair centre.
        double lead = 0.75 * kSidDay;
        double ahead = rel + rate * lead + 0.5 * A * lead * lead;
        if (sgn * rate > 0 && sgn * ahead >= half_box) {
            double dv = r0 / 3.0 * (rate + sgn * v_target);   // lam_dot changes by -3 dv / a
            half = dv / 2;
            burn(half);
            pending_half = per_day / 2;
            lockout = 3 * per_day;
            if (first_burn >= 0) dv_cycles += std::fabs(dv);
            else first_burn = t;
            last_burn = t;
            ++pairs;
        }
    }
    if (pairs > 1) {
        out.dv_rate = dv_cycles / (last_burn - first_burn) * kYear;
        out.cycle_days = (last_burn - first_burn) / (pairs - 1) / kDay;
    }
    return out;
}

}  // namespace geo
