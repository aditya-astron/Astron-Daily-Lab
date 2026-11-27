// Astron Built this AP
// Project 067 | Attitude control | Gravity-gradient stability regions and eccentric-orbit pitch libration
// Author: Aditya (aditya-astron)
// Date: 2026-11-27 IST
//
// Planar pitch in an eccentric orbit, true anomaly as the independent
// variable (' = d/dnu), pitch theta about the body y axis (y = -orbit normal):
//
//   (1 + e cos nu) theta'' - 2 e sin nu (theta' - 1) + 3 k2 sin theta cos theta = 0
//
// The "- 1" is the orbit frame's own non-uniform rotation: the body has no
// reason to follow the swing of the local vertical, so eccentricity forces
// pitch even at zero initial error. The sign of that term is specific to
// this frame and is pinned by comparison with the 3D model.
//
// Linearised and homogeneous, the equation is periodic in nu with zero mean
// damping (int 2e sin nu/(1+e cos nu) over 2pi = 0), so det(monodromy) = 1 and
// the motion is bounded iff |tr Phi(2pi)| < 2. Removing the first-derivative
// term with theta = u exp(-e cos nu) gives, to first order in e, a Mathieu
// equation u'' + (delta + eps cos nu) u = 0 with
//   delta = 3 k2,   eps = e (1 - 3 k2).
// Parametric tongues sit at sqrt(3 k2) = n/2, i.e. k2 = n^2/12; the first
// (n = 1) has width delta = 1/4 +- eps/2, so k2 in 1/12 -+ e/8. Higher
// tongues open as eps^n. The n = 2 tongue never opens: at k2 = 1/3 the
// monodromy is the identity for every e (both solutions 2pi-periodic), so
// that point is a pure forced resonance, not a parametric one.

use crate::ode::{rk4, transition_matrix};

pub const STEPS_PER_ORBIT: usize = 4000;

pub fn rhs_nonlinear(k2: f64, e: f64) -> impl Fn(f64, &[f64], &mut [f64]) {
    move |nu, x, dx| {
        let (s, c) = nu.sin_cos();
        dx[0] = x[1];
        dx[1] = (2.0 * e * s * (x[1] - 1.0) - 1.5 * k2 * (2.0 * x[0]).sin()) / (1.0 + e * c);
    }
}

/// Linear homogeneous pitch about theta_ref(nu) = 0 (or about a reference
/// trajectory when `cos2` is supplied by the caller, see `monodromy_about`).
pub fn rhs_linear(k2: f64, e: f64) -> impl Fn(f64, &[f64], &mut [f64]) {
    move |nu, x, dx| {
        let (s, c) = nu.sin_cos();
        dx[0] = x[1];
        dx[1] = (2.0 * e * s * x[1] - 3.0 * k2 * x[0]) / (1.0 + e * c);
    }
}

/// Monodromy matrix of the homogeneous linear pitch equation, row-major 2x2.
pub fn monodromy(k2: f64, e: f64, steps: usize) -> [f64; 4] {
    let p = transition_matrix(
        &rhs_linear(k2, e),
        2,
        0.0,
        2.0 * std::f64::consts::PI,
        steps,
    );
    [p[0], p[1], p[2], p[3]]
}

pub fn trace(k2: f64, e: f64) -> f64 {
    let m = monodromy(k2, e, STEPS_PER_ORBIT);
    m[0] + m[3]
}

/// Largest Floquet multiplier magnitude (det = 1 so the pair is lambda, 1/lambda
/// when |tr| > 2 and a unit-circle pair otherwise).
pub fn max_multiplier(tr: f64) -> f64 {
    if tr.abs() <= 2.0 {
        1.0
    } else {
        let a = tr.abs();
        0.5 * (a + (a * a - 4.0).sqrt())
    }
}

/// Edges of the n-th parametric tongue near k2 = n^2/12 at eccentricity e.
/// Inside, (-1)^n tr > 2. Returns None if the tongue is closed to within the
/// search resolution.
pub fn tongue(n: u32, e: f64) -> Option<(f64, f64)> {
    let centre = (n * n) as f64 / 12.0;
    let sgn = if n % 2 == 0 { 1.0 } else { -1.0 };
    let f = |k: f64| sgn * trace(k, e) - 2.0;
    let half = 0.02 + 0.5 * e;
    let (mut a, mut b) = ((centre - half).max(1e-6), centre + half);

    // golden-section for the peak of sgn*tr
    let g = 0.5 * (5f64.sqrt() - 1.0);
    let mut x1 = b - g * (b - a);
    let mut x2 = a + g * (b - a);
    let (mut f1, mut f2) = (f(x1), f(x2));
    for _ in 0..80 {
        if f1 > f2 {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = b - g * (b - a);
            f1 = f(x1);
        } else {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + g * (b - a);
            f2 = f(x2);
        }
        if b - a < 1e-13 {
            break;
        }
    }
    let peak = 0.5 * (a + b);
    // tr is good to ~1e-14 at 4000 steps, so an excess under 5e-14 is round-off.
    // Inside a tongue of width w the excess is only ~ (d2 tr/dk2) w^2 / 8, which
    // limits resolvable widths to a few 1e-8.
    if f(peak) <= 5e-14 {
        return None;
    }
    let bisect = |mut inside: f64, mut outside: f64| {
        for _ in 0..100 {
            let m = 0.5 * (inside + outside);
            if f(m) > 0.0 {
                inside = m;
            } else {
                outside = m;
            }
            if (inside - outside).abs() < 1e-14 {
                break;
            }
        }
        0.5 * (inside + outside)
    };
    let lo = bisect(peak, (centre - half).max(1e-6));
    let hi = bisect(peak, (centre + half).min(1.0));
    Some((lo, hi))
}

/// Flow of the nonlinear equation over one orbit plus its 2x2 variational
/// matrix, integrated together.
pub fn flow_with_jacobian(k2: f64, e: f64, x0: [f64; 2], steps: usize) -> ([f64; 2], [f64; 4]) {
    let f = move |nu: f64, y: &[f64], dy: &mut [f64]| {
        let (s, c) = nu.sin_cos();
        let d = 1.0 + e * c;
        dy[0] = y[1];
        dy[1] = (2.0 * e * s * (y[1] - 1.0) - 1.5 * k2 * (2.0 * y[0]).sin()) / d;
        // variational: d'' = (2 e s d' - 3 k2 cos(2 theta) d) / (1 + e c)
        let kk = 3.0 * k2 * (2.0 * y[0]).cos();
        for col in 0..2 {
            let (a, b) = (y[2 + col], y[4 + col]);
            dy[2 + col] = b;
            dy[4 + col] = (2.0 * e * s * b - kk * a) / d;
        }
    };
    let mut y = vec![x0[0], x0[1], 1.0, 0.0, 0.0, 1.0];
    rk4(
        &f,
        0.0,
        2.0 * std::f64::consts::PI,
        &mut y,
        steps,
        |_, _| {},
    );
    ([y[0], y[1]], [y[2], y[3], y[4], y[5]])
}

/// 2pi-periodic pitch libration by Newton shooting, continued from the
/// linear forced solution. Returns (theta0, theta0', monodromy).
pub fn periodic_orbit(k2: f64, e: f64, guess: [f64; 2]) -> Option<([f64; 2], [f64; 4])> {
    let mut x = guess;
    for _ in 0..50 {
        let (xf, m) = flow_with_jacobian(k2, e, x, STEPS_PER_ORBIT);
        let r = [xf[0] - x[0], xf[1] - x[1]];
        if r[0].abs() + r[1].abs() < 1e-12 {
            return Some((x, m));
        }
        let j = [m[0] - 1.0, m[1], m[2], m[3] - 1.0];
        let det = j[0] * j[3] - j[1] * j[2];
        if det.abs() < 1e-14 {
            return None;
        }
        let dx0 = (j[3] * r[0] - j[1] * r[1]) / det;
        let dx1 = (-j[2] * r[0] + j[0] * r[1]) / det;
        x = [x[0] - dx0, x[1] - dx1];
        if !x[0].is_finite() || x[0].abs() > 3.0 {
            return None;
        }
    }
    None
}

/// First-order forced amplitude: theta ~ A sin(nu) with A = -2e/(3k2 - 1).
pub fn forced_amplitude_linear(k2: f64, e: f64) -> f64 {
    -2.0 * e / (3.0 * k2 - 1.0)
}

/// Peak |theta| over one orbit of a trajectory starting at x0.
pub fn peak_over_orbit(k2: f64, e: f64, x0: [f64; 2]) -> f64 {
    let mut y = vec![x0[0], x0[1]];
    let mut pk = x0[0].abs();
    rk4(
        &rhs_nonlinear(k2, e),
        0.0,
        2.0 * std::f64::consts::PI,
        &mut y,
        STEPS_PER_ORBIT,
        |_, s| pk = pk.max(s[0].abs()),
    );
    pk
}
