// Astron Built this AP
// Project 047 | Navigation | GNSS snapshot positioning with RAIM fault detection and exclusion
// Author: Aditya (aditya-astron)
// Date: 2026-11-19 IST
//
// Study for a receiver in Mumbai on a 24-satellite GPS-like constellation:
//   1. one-epoch fix, DOPs and per-satellite slopes / protection levels
//   2. a ramped bias on the hardest-to-see satellite: detection and exclusion
//   3. 24 h RAIM and FDE availability against a 556 m HAL (0.3 nmi, NPA),
//      for a modern 5 m and an SA-era 33.3 m pseudorange sigma, with the
//      full constellation and with each single satellite out.

#include <algorithm>
#include <cstdio>
#include <random>

#include "raim.hpp"

using namespace raim;

namespace {

const double kLat = rad(19.0760), kLon = rad(72.8777), kH = 14.0;
const double kMask = rad(5.0);

struct Sky {
    std::vector<int> ids;
    std::vector<Vec3> pos;
    std::vector<double> az, el;
};

Sky visible(const std::vector<Slot>& c, double t, const Vec3& rx, int skip = -1) {
    Sky s;
    for (int i = 0; i < static_cast<int>(c.size()); ++i) {
        if (i == skip) continue;
        const Vec3 p = sat_ecef(c[i], t);
        double az, el;
        az_el(rx, kLat, kLon, p, az, el);
        if (el < kMask) continue;
        s.ids.push_back(i);
        s.pos.push_back(p);
        s.az.push_back(az);
        s.el.push_back(el);
    }
    return s;
}

Mat sky_h(const Sky& s) { return geometry_enu(s.az, s.el); }

// Detection-only HPL of every N-1 subset; FDE needs all of them under HAL.
bool fde_available(const Sky& s, const RaimConfig& cfg, double hal) {
    if (s.ids.size() < 6) return false;
    const Mat H = sky_h(s);
    for (size_t i = 0; i < H.size(); ++i) {
        Mat Hs;
        for (size_t j = 0; j < H.size(); ++j)
            if (j != i) Hs.push_back(H[j]);
        if (raim_geometry(Hs, cfg).hpl > hal) return false;
    }
    return true;
}

const char* plane_slot(int id) {
    static char buf[8];
    std::snprintf(buf, sizeof buf, "%c%d", 'A' + id / 4, id % 4 + 1);
    return buf;
}

}  // namespace

int main() {
    const auto constel = gps_like_24();
    const Vec3 rx = geodetic_to_ecef(kLat, kLon, kH);
    const double clock_bias = 85e3;  // m, about 0.28 ms of receiver clock offset
    const double hal = 556.0;

    // ---------------------------------------------------------- 1. snapshot
    const double t0 = 0.0;
    const Sky sky = visible(constel, t0, rx);
    std::vector<double> rho;
    for (const auto& p : sky.pos) rho.push_back(true_pseudorange(p, rx, clock_bias));
    const Fix fix = solve(sky.pos, rho);
    const Mat H = to_enu(fix.H, kLat, kLon);
    const RaimConfig cfg5{5.0, 1e-5, 1e-3};
    const RaimGeometry g = raim_geometry(H, cfg5);
    const Dop d = dop(H);

    std::printf("1. Snapshot, Mumbai 19.076 N 72.878 E, t = 0, mask 5 deg\n");
    std::printf("   %zu satellites in view; noiseless fix from Earth's centre in %d iterations,\n",
                sky.ids.size(), fix.iterations);
    std::printf("   position error %.2e m, clock error %.2e m\n", norm(sub(fix.r, rx)), fix.b - clock_bias);
    std::printf("   GDOP %.2f  PDOP %.2f  HDOP %.2f  VDOP %.2f  TDOP %.2f\n", d.gdop, d.pdop, d.hdop, d.vdop, d.tdop);
    std::printf("   sigma 5 m, Pfa 1e-5, Pmd 1e-3: T = %.3f (chi2, %d dof), lambda_det = %.3f\n", g.threshold, g.dof,
                g.lambda_det);
    std::printf("   %-5s %7s %6s %7s %7s %7s\n", "sat", "az", "el", "S_ii", "Hslope", "Vslope");
    for (size_t i = 0; i < sky.ids.size(); ++i)
        std::printf("   %-5s %7.1f %6.1f %7.3f %7.3f %7.3f%s\n", plane_slot(sky.ids[i]), deg(sky.az[i]),
                    deg(sky.el[i]), g.s_diag[i], g.hslope[i], g.vslope[i],
                    static_cast<int>(i) == g.worst_h ? "  <- worst" : "");
    std::printf("   HPL %.1f m  VPL %.1f m  (sigma 33.3 m: HPL %.1f m)\n\n", g.hpl, g.vpl,
                raim_geometry(H, {33.3, 1e-5, 1e-3}).hpl);

    // ------------------------------------------------- 2. ramped fault demo
    std::printf("2. Bias ramp on %s (worst slope), sigma 5 m, one noise draw per step\n", plane_slot(sky.ids[g.worst_h]));
    std::printf("   %8s %9s %8s %6s %9s %12s\n", "bias m", "SSE/s^2", "horiz m", "alarm", "excluded", "horiz after m");
    std::mt19937_64 rng(47);
    std::normal_distribution<double> nd(0.0, cfg5.sigma);
    const auto B = enu_basis(kLat, kLon);
    auto horiz = [&](const Vec3& r) {
        const Vec3 e = sub(r, rx);
        const double de = e[0] * B[0][0] + e[1] * B[0][1] + e[2] * B[0][2];
        const double dn = e[0] * B[1][0] + e[1] * B[1][1] + e[2] * B[1][2];
        return std::hypot(de, dn);
    };
    for (double bias : {0.0, 25.0, 50.0, 75.0, 100.0, 150.0, 250.0}) {
        std::vector<double> r = rho;
        for (double& x : r) x += nd(rng);
        r[g.worst_h] += bias;
        const FdeResult fr = detect_and_exclude(sky.pos, r, cfg5);
        const double h0 = horiz(solve(sky.pos, r).r);
        double h1 = h0;
        if (fr.recovered) {
            std::vector<Vec3> s2;
            std::vector<double> r2;
            for (size_t j = 0; j < r.size(); ++j)
                if (static_cast<int>(j) != fr.excluded) {
                    s2.push_back(sky.pos[j]);
                    r2.push_back(r[j]);
                }
            h1 = horiz(solve(s2, r2).r);
        }
        std::printf("   %8.0f %9.2f %8.2f %6s %9s %12.2f\n", bias, fr.test, h0, fr.alarm ? "yes" : "no",
                    fr.recovered ? plane_slot(sky.ids[fr.excluded]) : "-", h1);
    }
    std::printf("   bias at which the test misses 1e-3 of the time: %.1f m -> mean horizontal error %.1f m = HPL\n\n",
                cfg5.sigma * std::sqrt(g.lambda_det / g.s_diag[g.worst_h]), g.hpl);

    // ----------------------------------------------- 3. 24 h availability
    std::printf("3. 24 h availability, 60 s steps, HAL %.0f m, Pfa 1e-5, Pmd 1e-3\n", hal);
    const int steps = 1440;
    for (double sigma : {5.0, 33.3}) {
        const RaimConfig cfg{sigma, 1e-5, 1e-3};
        // Full constellation, then every single-satellite outage.
        int worst_out = -1, worst_det = steps + 1, worst_fde = steps + 1;
        int full_det = 0, full_fde = 0, nmin = 99, nmax = 0;
        double pdop_max = 0.0;
        std::vector<double> hpls;
        for (int skip = -1; skip < 24; ++skip) {
            int det = 0, fde = 0;
            for (int k = 0; k < steps; ++k) {
                const Sky s = visible(constel, 60.0 * k, rx, skip);
                if (s.ids.size() < 5) continue;
                const Mat Hs = sky_h(s);
                const RaimGeometry gs = raim_geometry(Hs, cfg);
                if (gs.hpl <= hal) ++det;
                if (fde_available(s, cfg, hal)) ++fde;
                if (skip == -1) {
                    nmin = std::min<int>(nmin, s.ids.size());
                    nmax = std::max<int>(nmax, s.ids.size());
                    pdop_max = std::max(pdop_max, dop(Hs).pdop);
                    hpls.push_back(gs.hpl);
                }
            }
            if (skip == -1) {
                full_det = det;
                full_fde = fde;
            } else if (det < worst_det || (det == worst_det && fde < worst_fde)) {
                worst_det = det;
                worst_fde = fde;
                worst_out = skip;
            }
        }
        std::sort(hpls.begin(), hpls.end());
        std::printf("   sigma %4.1f m: %d-%d in view, max PDOP %.2f, HPL median %.0f m, 99th %.0f m, max %.0f m\n",
                    sigma, nmin, nmax, pdop_max, hpls[hpls.size() / 2], hpls[hpls.size() * 99 / 100], hpls.back());
        std::printf("      all 24 up:        detection %6.2f %%   FDE %6.2f %%\n", 100.0 * full_det / steps,
                    100.0 * full_fde / steps);
        std::printf("      worst one out %s: detection %6.2f %%   FDE %6.2f %%  (%d / %d min lost)\n",
                    plane_slot(worst_out), 100.0 * worst_det / steps, 100.0 * worst_fde / steps, steps - worst_det,
                    steps - worst_fde);
    }
    return 0;
}
