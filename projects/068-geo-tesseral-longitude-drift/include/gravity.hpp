// Astron Built this AP
// Project 068 | Numerical methods / Orbital mechanics | Spherical-harmonic gravity and GEO longitude drift
// Author: Aditya (aditya-astron)
// Date: 2026-11-28 IST
//
// Spherical-harmonic gravity field, Earth-fixed frame, fully normalized
// (4 pi) coefficients without the Condon-Shortley phase:
//
//   U = GM/r sum_n (R/r)^n sum_m Pbar_nm(sin phi) (C_nm cos m lam + S_nm sin m lam)
//
// Two independent acceleration paths so each can check the other:
//   accel_spherical  normalized Legendre functions + spherical partials
//                    (singular only at the poles, stable to high degree)
//   accel_cunningham Montenbruck & Gill V/W recursion on unnormalized
//                    coefficients (non-singular, overflows past n ~ 150)

#pragma once
#include <array>
#include <cmath>
#include <vector>

namespace grav {

using Vec3 = std::array<double, 3>;

struct Field {
    int N = 0;
    double GM = 0.0, R = 0.0;
    std::vector<double> C, S;   // packed lower triangle, index tri(n, m)

    Field(int nmax, double gm, double radius)
        : N(nmax), GM(gm), R(radius),
          C(tri(nmax + 1, 0), 0.0), S(tri(nmax + 1, 0), 0.0) { C[0] = 1.0; }

    static int tri(int n, int m) { return n * (n + 1) / 2 + m; }
    double& c(int n, int m) { return C[tri(n, m)]; }
    double& s(int n, int m) { return S[tri(n, m)]; }
    double c(int n, int m) const { return C[tri(n, m)]; }
    double s(int n, int m) const { return S[tri(n, m)]; }

    // Truncated copy, used to isolate the effect of each degree.
    Field truncated(int nmax) const {
        Field f(nmax, GM, R);
        for (int n = 0; n <= nmax; ++n)
            for (int m = 0; m <= n; ++m) { f.c(n, m) = c(n, m); f.s(n, m) = s(n, m); }
        return f;
    }
};

// EGM96 coefficients through degree and order 4 (fully normalized).
inline Field egm96_4x4() {
    Field f(4, 3.986004415e14, 6378136.3);
    f.c(2, 0) = -4.84165371736e-04;
    f.c(2, 1) = -1.86987635955e-10; f.s(2, 1) =  1.19528012031e-09;
    f.c(2, 2) =  2.43914352398e-06; f.s(2, 2) = -1.40016683654e-06;
    f.c(3, 0) =  9.57254173792e-07;
    f.c(3, 1) =  2.03046201047e-06; f.s(3, 1) =  2.48200415856e-07;
    f.c(3, 2) =  9.04787894809e-07; f.s(3, 2) = -6.19005475177e-07;
    f.c(3, 3) =  7.21321757121e-07; f.s(3, 3) =  1.41434926192e-06;
    f.c(4, 0) =  5.39873863789e-07;
    f.c(4, 1) = -5.36321616971e-07; f.s(4, 1) = -4.73440265853e-07;
    f.c(4, 2) =  3.50694105785e-07; f.s(4, 2) =  6.62671572540e-07;
    f.c(4, 3) =  9.90771803829e-07; f.s(4, 3) = -2.00928369177e-07;
    f.c(4, 4) = -1.88560802735e-07; f.s(4, 4) =  3.08853169333e-07;
    return f;
}

// Normalization Pbar_nm = N_nm P_nm,  N_nm = sqrt((2 - d_m0)(2n+1)(n-m)!/(n+m)!)
inline double norm_factor(int n, int m) {
    double lr = std::lgamma(n - m + 1.0) - std::lgamma(n + m + 1.0);
    return std::sqrt((m == 0 ? 1.0 : 2.0) * (2 * n + 1) * std::exp(lr));
}

// Fully normalized Pbar_nm(t), t = sin(phi), u = cos(phi), for n <= nmax.
// Sectoral seed then the standard two-term column recursion in n; every
// coefficient is O(1), so it stays accurate far past where P_nm overflows.
inline void legendre(int nmax, double t, double u, std::vector<double>& P) {
    P.assign(Field::tri(nmax + 1, 0), 0.0);
    P[0] = 1.0;
    if (nmax == 0) return;
    P[Field::tri(1, 1)] = std::sqrt(3.0) * u;
    for (int m = 2; m <= nmax; ++m)
        P[Field::tri(m, m)] = u * std::sqrt((2.0 * m + 1.0) / (2.0 * m)) * P[Field::tri(m - 1, m - 1)];
    for (int m = 0; m <= nmax; ++m) {
        if (m + 1 > nmax) break;
        P[Field::tri(m + 1, m)] = std::sqrt(2.0 * m + 3.0) * t * P[Field::tri(m, m)];
        for (int n = m + 2; n <= nmax; ++n) {
            double a = std::sqrt((2.0 * n - 1) * (2.0 * n + 1) / double((n - m) * (n + m)));
            double b = std::sqrt((2.0 * n + 1) * (n + m - 1) * (n - m - 1) /
                                 double((n - m) * (n + m) * (2.0 * n - 3)));
            P[Field::tri(n, m)] = a * t * P[Field::tri(n - 1, m)] - b * P[Field::tri(n - 2, m)];
        }
    }
}

struct Spherical { double r, phi, lam; };

inline Spherical to_spherical(const Vec3& x) {
    double rho = std::hypot(x[0], x[1]);
    return {std::hypot(rho, x[2]), std::atan2(x[2], rho), std::atan2(x[1], x[0])};
}

inline double potential(const Field& f, const Vec3& x) {
    auto [r, phi, lam] = to_spherical(x);
    std::vector<double> P;
    legendre(f.N, std::sin(phi), std::cos(phi), P);
    double sum = 0.0, q = 1.0, ratio = f.R / r;
    for (int n = 0; n <= f.N; ++n, q *= ratio) {
        double inner = 0.0;
        for (int m = 0; m <= n; ++m)
            inner += P[Field::tri(n, m)] * (f.c(n, m) * std::cos(m * lam) + f.s(n, m) * std::sin(m * lam));
        sum += q * inner;
    }
    return f.GM / r * sum;
}

// Gradient of U from the spherical partials. Uses
//   dPbar_nm/dphi = k_nm Pbar_n,m+1 - m tan(phi) Pbar_nm,
//   k_nm = sqrt(n(n+1)/2) for m = 0, sqrt((n-m)(n+m+1)) otherwise,
// which follows from dP_nm/dphi = P_n,m+1 - m tan(phi) P_nm and N_nm/N_n,m+1.
// Returns the local (radial, north, east) components.
inline Vec3 accel_local(const Field& f, double r, double phi, double lam) {
    double t = std::sin(phi), u = std::cos(phi), tanp = t / u;
    std::vector<double> P;
    legendre(f.N, t, u, P);

    double dUr = 0.0, dUphi = 0.0, dUlam = 0.0, q = 1.0, ratio = f.R / r;
    for (int n = 0; n <= f.N; ++n, q *= ratio) {
        double sr = 0.0, sp = 0.0, sl = 0.0;
        for (int m = 0; m <= n; ++m) {
            double cm = std::cos(m * lam), sm = std::sin(m * lam);
            double Cnm = f.c(n, m), Snm = f.s(n, m);
            double Pnm = P[Field::tri(n, m)];
            double Pnm1 = (m < n) ? P[Field::tri(n, m + 1)] : 0.0;
            double k = (m == 0) ? std::sqrt(n * (n + 1) / 2.0) : std::sqrt(double(n - m) * (n + m + 1));
            double dP = k * Pnm1 - m * tanp * Pnm;
            double cs = Cnm * cm + Snm * sm;
            sr += Pnm * cs;
            sp += dP * cs;
            sl += m * Pnm * (Snm * cm - Cnm * sm);
        }
        dUr -= (n + 1) * q * sr;
        dUphi += q * sp;
        dUlam += q * sl;
    }
    double g = f.GM / r;
    dUr *= g / r;
    dUphi *= g;
    dUlam *= g;
    return {dUr, dUphi / r, dUlam / (r * u)};
}

inline Vec3 accel_spherical(const Field& f, const Vec3& x) {
    auto [r, phi, lam] = to_spherical(x);
    Vec3 loc = accel_local(f, r, phi, lam);
    double t = std::sin(phi), u = std::cos(phi);
    double cl = std::cos(lam), sl = std::sin(lam);
    Vec3 er{u * cl, u * sl, t}, ep{-t * cl, -t * sl, u}, el{-sl, cl, 0.0};
    Vec3 a{};
    for (int i = 0; i < 3; ++i)
        a[i] = loc[0] * er[i] + loc[1] * ep[i] + loc[2] * el[i];
    return a;
}

// Montenbruck & Gill, Satellite Orbits, eqs. 3.29-3.33.
inline Vec3 accel_cunningham(const Field& f, const Vec3& x) {
    const int M = f.N + 1;
    const int sz = Field::tri(M + 1, 0);
    std::vector<double> V(sz, 0.0), W(sz, 0.0);
    double r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
    double x0 = f.R * x[0] / r2, y0 = f.R * x[1] / r2, z0 = f.R * x[2] / r2, rho = f.R * f.R / r2;
    auto id = Field::tri;

    V[0] = f.R / std::sqrt(r2);
    W[0] = 0.0;
    for (int m = 0; m <= M; ++m) {
        if (m > 0) {
            V[id(m, m)] = (2 * m - 1) * (x0 * V[id(m - 1, m - 1)] - y0 * W[id(m - 1, m - 1)]);
            W[id(m, m)] = (2 * m - 1) * (x0 * W[id(m - 1, m - 1)] + y0 * V[id(m - 1, m - 1)]);
        }
        if (m + 1 <= M) {
            V[id(m + 1, m)] = (2 * m + 1) * z0 * V[id(m, m)];
            W[id(m + 1, m)] = (2 * m + 1) * z0 * W[id(m, m)];
        }
        for (int n = m + 2; n <= M; ++n) {
            double a = double(2 * n - 1) / (n - m), b = double(n + m - 1) / (n - m);
            V[id(n, m)] = a * z0 * V[id(n - 1, m)] - b * rho * V[id(n - 2, m)];
            W[id(n, m)] = a * z0 * W[id(n - 1, m)] - b * rho * W[id(n - 2, m)];
        }
    }

    double ax = 0.0, ay = 0.0, az = 0.0;
    for (int n = 0; n <= f.N; ++n) {
        for (int m = 0; m <= n; ++m) {
            double Nf = norm_factor(n, m);
            double C = f.c(n, m) * Nf, S = f.s(n, m) * Nf;
            if (m == 0) {
                ax -= C * V[id(n + 1, 1)];
                ay -= C * W[id(n + 1, 1)];
            } else {
                double fac = double(n - m + 2) * (n - m + 1);
                ax += 0.5 * (-C * V[id(n + 1, m + 1)] - S * W[id(n + 1, m + 1)]
                             + fac * (C * V[id(n + 1, m - 1)] + S * W[id(n + 1, m - 1)]));
                ay += 0.5 * (-C * W[id(n + 1, m + 1)] + S * V[id(n + 1, m + 1)]
                             + fac * (-C * W[id(n + 1, m - 1)] + S * V[id(n + 1, m - 1)]));
            }
            az += (n - m + 1) * (-C * V[id(n + 1, m)] - S * W[id(n + 1, m)]);
        }
    }
    double g = f.GM / (f.R * f.R);
    return {g * ax, g * ay, g * az};
}

}  // namespace grav
