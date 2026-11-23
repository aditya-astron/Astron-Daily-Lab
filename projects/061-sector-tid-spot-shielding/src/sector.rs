// Astron Built this AP
// Project 061 | Radiation shielding | Sector-based TID analysis and spot-shield sizing for an avionics box
// Author: Aditya (aditya-astron)
// Date: 2026-11-23 IST
//
// Sector analysis: the environment is taken as isotropic, so the dose at a
// point is the solid-angle average of the solid-sphere dose-depth curve
// evaluated at the areal density seen along each direction,
//
//     D(p) = (1/4pi) * integral D_sphere(t(p, u)) dOmega(u).
//
// The sphere of directions is sampled with a Fibonacci lattice: N points of
// equal area, z_i = 1 - (2i+1)/N, phi_i = i * golden angle. In z it is the
// midpoint rule, so axisymmetric smooth integrands converge as 1/N^2; the
// edge and corner kinks of a box model cut that to roughly 1/N.

use crate::dose::DoseDepth;
use crate::geometry::{Scene, V3};

pub fn fibonacci_dirs(n: usize) -> Vec<V3> {
    let ga = std::f64::consts::PI * (3.0 - 5.0_f64.sqrt());
    (0..n)
        .map(|i| {
            let z = 1.0 - (2.0 * i as f64 + 1.0) / n as f64;
            let r = (1.0 - z * z).max(0.0).sqrt();
            let phi = ga * i as f64;
            [r * phi.cos(), r * phi.sin(), z]
        })
        .collect()
}

/// Uniform random directions from a fixed-seed xorshift64*, for the Monte
/// Carlo comparison.
pub fn random_dirs(n: usize, seed: u64) -> Vec<V3> {
    let mut s = seed.max(1);
    let mut next = || {
        s ^= s >> 12;
        s ^= s << 25;
        s ^= s >> 27;
        (s.wrapping_mul(0x2545_F491_4F6C_DD1D) >> 11) as f64 / (1u64 << 53) as f64
    };
    (0..n)
        .map(|_| {
            let z = 2.0 * next() - 1.0;
            let phi = 2.0 * std::f64::consts::PI * next();
            let r = (1.0 - z * z).sqrt();
            [r * phi.cos(), r * phi.sin(), z]
        })
        .collect()
}

pub fn thickness_map(scene: &Scene, p: &V3, dirs: &[V3]) -> Vec<f64> {
    dirs.iter().map(|u| scene.areal_density(p, u)).collect()
}

pub fn fold(dd: &DoseDepth, t: &[f64]) -> f64 {
    t.iter().map(|&x| dd.eval(x)).sum::<f64>() / t.len() as f64
}

pub fn sector_dose(scene: &Scene, p: &V3, dd: &DoseDepth, dirs: &[V3]) -> f64 {
    fold(dd, &thickness_map(scene, p, dirs))
}

/// Adaptive Simpson on [a, b] with the usual Richardson correction.
pub fn adaptive_simpson<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64, tol: f64) -> f64 {
    fn rec<F: Fn(f64) -> f64>(
        f: &F, a: f64, b: f64, fa: f64, fm: f64, fb: f64, whole: f64, tol: f64, depth: u32,
    ) -> f64 {
        let m = 0.5 * (a + b);
        let lm = 0.5 * (a + m);
        let rm = 0.5 * (m + b);
        let flm = f(lm);
        let frm = f(rm);
        let left = (m - a) / 6.0 * (fa + 4.0 * flm + fm);
        let right = (b - m) / 6.0 * (fm + 4.0 * frm + fb);
        let delta = left + right - whole;
        if depth == 0 || delta.abs() <= 15.0 * tol {
            return left + right + delta / 15.0;
        }
        rec(f, a, m, fa, flm, fm, left, 0.5 * tol, depth - 1)
            + rec(f, m, b, fm, frm, fb, right, 0.5 * tol, depth - 1)
    }
    let fa = f(a);
    let fb = f(b);
    let fm = f(0.5 * (a + b));
    let whole = (b - a) / 6.0 * (fa + 4.0 * fm + fb);
    rec(f, a, b, fa, fm, fb, whole, tol, 40)
}

/// Dose at the centre of a cubic shell of wall areal density w, by direct
/// quadrature over one eighth of a face. With the face at unit distance and
/// (x, y) on it, r = sqrt(1 + x^2 + y^2): the ray crosses the wall with
/// chord w*r and dOmega = dx dy / r^3. 6 faces * 4 quarters / (4 pi).
pub fn cube_centre_reference(dd: &DoseDepth, w: f64, tol: f64) -> f64 {
    let inner = |x: f64| {
        let g = |y: f64| {
            let r = (1.0 + x * x + y * y).sqrt();
            dd.eval(w * r) / (r * r * r)
        };
        adaptive_simpson(&g, 0.0, 1.0, tol)
    };
    24.0 / (4.0 * std::f64::consts::PI) * adaptive_simpson(&inner, 0.0, 1.0, tol)
}

/// Dose at distance d from the centre of a spherical shell (R1 < R2, d < R1,
/// density rho), by 1D quadrature in mu = cos(angle to the radius vector):
/// t(mu) = rho * (sqrt(R2^2 - d^2 + d^2 mu^2) - sqrt(R1^2 - d^2 + d^2 mu^2)).
pub fn offset_shell_reference(dd: &DoseDepth, r1: f64, r2: f64, rho: f64, d: f64, tol: f64) -> f64 {
    let f = |mu: f64| {
        let s = d * d * (mu * mu - 1.0);
        let t = rho * ((r2 * r2 + s).sqrt() - (r1 * r1 + s).sqrt());
        dd.eval(t)
    };
    0.5 * adaptive_simpson(&f, -1.0, 1.0, tol)
}

/// Linear-interpolated percentile of a sample (q in [0, 1]).
pub fn percentile(t: &[f64], q: f64) -> f64 {
    let mut v = t.to_vec();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let x = q * (v.len() - 1) as f64;
    let i = x.floor() as usize;
    let j = (i + 1).min(v.len() - 1);
    v[i] + (x - i as f64) * (v[j] - v[i])
}

/// Bisection for the root of a function that is decreasing in s on [lo, hi]:
/// returns s with f(s) = target, or None if the bracket does not hold it.
pub fn bisect_decreasing<F: Fn(f64) -> f64>(f: F, target: f64, mut lo: f64, mut hi: f64, tol: f64) -> Option<f64> {
    if f(lo) <= target {
        return Some(lo);
    }
    if f(hi) > target {
        return None;
    }
    while hi - lo > tol {
        let mid = 0.5 * (lo + hi);
        if f(mid) > target { lo = mid } else { hi = mid }
    }
    Some(0.5 * (lo + hi))
}
