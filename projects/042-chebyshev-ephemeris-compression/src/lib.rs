// Astron Built this AP
// Project 042 | Numerical methods / Navigation | Chebyshev ephemeris compression for an eccentric orbit
// Author: Aditya (aditya-astron)
// Date: 2026-11-16 IST

//! Piecewise Chebyshev ephemeris in the style of the JPL DE files and the
//! CCSDS/NAIF type-2 SPK segments: each time segment stores N+1 Chebyshev
//! coefficients per position axis, velocity comes from differentiating the
//! same polynomial, and the segment layout is chosen to hit a position
//! tolerance with as few stored doubles as possible.
//!
//! The truth model is an exact two-body Kepler orbit, so every error printed
//! is interpolation error and nothing else. That makes the convergence rate
//! predictable: x(t) is analytic except where dM/dE = 1 - e cos E = 0, i.e. at
//! E = 2 pi m +/- i acosh(1/e). Those branch points sit in the complex
//! mean-anomaly plane at M* = 2 pi m +/- i eta with
//!     eta = acosh(1/e) - sqrt(1 - e^2),
//! and Chebyshev coefficients on a segment decay like rho^-k, where rho is the
//! Bernstein ellipse through the nearest M* after mapping the segment to [-1, 1].

pub mod complex;
use complex::C64;
use std::f64::consts::PI;

pub const MU_EARTH: f64 = 398_600.4418; // km^3/s^2

// ---------------------------------------------------------------- Chebyshev

/// Chebyshev points of the first kind on [-1, 1], ordered from +1 to -1.
pub fn cheb_nodes(n: usize) -> Vec<f64> {
    let m = (n + 1) as f64;
    (0..=n).map(|j| (PI * (j as f64 + 0.5) / m).cos()).collect()
}

/// Coefficients of the degree-n interpolant at first-kind nodes, using the
/// discrete orthogonality of T_k on those nodes (a type-II DCT, done
/// directly since n is small). Convention: f(x) = sum_{k=0}^{n} c_k T_k(x).
pub fn cheb_fit<F: Fn(f64) -> f64>(f: F, n: usize) -> Vec<f64> {
    let m = (n + 1) as f64;
    let fx: Vec<f64> = cheb_nodes(n).into_iter().map(&f).collect();
    let mut c: Vec<f64> = (0..=n)
        .map(|k| {
            let s: f64 = fx
                .iter()
                .enumerate()
                .map(|(j, &v)| v * (PI * k as f64 * (j as f64 + 0.5) / m).cos())
                .sum();
            2.0 * s / m
        })
        .collect();
    c[0] *= 0.5;
    c
}

/// Clenshaw recurrence: sum c_k T_k(x) in O(n) with good rounding behaviour.
pub fn cheb_eval(c: &[f64], x: f64) -> f64 {
    let (mut b1, mut b2) = (0.0, 0.0);
    for &ck in c.iter().skip(1).rev() {
        let b0 = 2.0 * x * b1 - b2 + ck;
        b2 = b1;
        b1 = b0;
    }
    x * b1 - b2 + c[0]
}

/// Coefficients of d/dx sum c_k T_k: d_{k-1} = d_{k+1} + 2 k c_k, halve d_0.
pub fn cheb_deriv(c: &[f64]) -> Vec<f64> {
    let n = c.len() - 1;
    if n == 0 {
        return vec![0.0];
    }
    let mut d = vec![0.0; n + 2];
    for k in (1..=n).rev() {
        d[k - 1] = d[k + 1] + 2.0 * k as f64 * c[k];
    }
    d[0] *= 0.5;
    d.truncate(n);
    d
}

// ------------------------------------------------------------------- Kepler

#[derive(Clone, Copy, Debug)]
pub struct Orbit {
    pub a: f64,    // km
    pub e: f64,    //
    pub inc: f64,  // rad
    pub raan: f64, // rad
    pub argp: f64, // rad
    pub mu: f64,   // km^3/s^2
}

impl Orbit {
    pub fn mean_motion(&self) -> f64 {
        (self.mu / self.a.powi(3)).sqrt()
    }
    pub fn period(&self) -> f64 {
        2.0 * PI / self.mean_motion()
    }

    /// Imaginary offset of the Kepler-equation branch points from the real
    /// mean-anomaly axis.
    pub fn eta(&self) -> f64 {
        (1.0 / self.e).acosh() - (1.0 - self.e * self.e).sqrt()
    }

    /// Inertial position (km) and velocity (km/s) at t seconds after perigee.
    pub fn state(&self, t: f64) -> ([f64; 3], [f64; 3]) {
        let n = self.mean_motion();
        let ea = solve_kepler(n * t, self.e);
        let (se, ce) = ea.sin_cos();
        let b = self.a * (1.0 - self.e * self.e).sqrt();
        let edot = n / (1.0 - self.e * ce);
        let p = [self.a * (ce - self.e), b * se];
        let v = [-self.a * se * edot, b * ce * edot];
        (self.rotate(p), self.rotate(v))
    }

    // Perifocal -> inertial: R3(-raan) R1(-inc) R3(-argp).
    fn rotate(&self, q: [f64; 2]) -> [f64; 3] {
        let (so, co) = self.raan.sin_cos();
        let (si, ci) = self.inc.sin_cos();
        let (sw, cw) = self.argp.sin_cos();
        let px = [co * cw - so * sw * ci, so * cw + co * sw * ci, sw * si];
        let py = [-co * sw - so * cw * ci, -so * sw + co * cw * ci, cw * si];
        [
            px[0] * q[0] + py[0] * q[1],
            px[1] * q[0] + py[1] * q[1],
            px[2] * q[0] + py[2] * q[1],
        ]
    }
}

/// Newton on E - e sin E = M with the Danby starting guess; converges to
/// machine precision in a handful of steps for e < 1.
pub fn solve_kepler(m: f64, e: f64) -> f64 {
    let m = m.rem_euclid(2.0 * PI);
    let mut ea = m + 0.85 * e * m.sin().signum();
    for _ in 0..50 {
        let f = ea - e * ea.sin() - m;
        let d = f / (1.0 - e * ea.cos());
        ea -= d;
        if d.abs() < 1e-15 {
            break;
        }
    }
    ea
}

// ----------------------------------------------------- Bernstein ellipse

/// Bernstein parameter rho of the ellipse with foci +/-1 passing through z.
pub fn bernstein_rho(z: C64) -> f64 {
    let w = (z - C64::new(1.0, 0.0)).sqrt() * (z + C64::new(1.0, 0.0)).sqrt();
    (z + w).abs().max((z - w).abs())
}

/// Predicted coefficient decay ratio for a segment [t0, t1] (s from perigee):
/// the smallest rho over the branch points M* = 2 pi m +/- i eta.
pub fn segment_rho(orb: &Orbit, t0: f64, t1: f64) -> f64 {
    let n = orb.mean_motion();
    let (mc, hm) = (0.5 * n * (t0 + t1), 0.5 * n * (t1 - t0));
    let eta = orb.eta();
    let k0 = (mc / (2.0 * PI)).round() as i64;
    let mut best = f64::INFINITY;
    for m in (k0 - 2)..=(k0 + 2) {
        for s in [1.0, -1.0] {
            let z = C64::new((2.0 * PI * m as f64 - mc) / hm, s * eta / hm);
            best = best.min(bernstein_rho(z));
        }
    }
    best
}

// ------------------------------------------------------------- Segments

/// One ephemeris record: time span and N+1 coefficients per axis.
#[derive(Clone, Debug)]
pub struct Segment {
    pub t0: f64,
    pub t1: f64,
    pub coef: [Vec<f64>; 3],
}

impl Segment {
    pub fn fit(orb: &Orbit, t0: f64, t1: f64, n: usize) -> Self {
        let (mid, half) = (0.5 * (t0 + t1), 0.5 * (t1 - t0));
        let axis = |i: usize| cheb_fit(|x| orb.state(mid + half * x).0[i], n);
        Segment { t0, t1, coef: [axis(0), axis(1), axis(2)] }
    }

    fn tau(&self, t: f64) -> f64 {
        (2.0 * t - self.t0 - self.t1) / (self.t1 - self.t0)
    }

    pub fn position(&self, t: f64) -> [f64; 3] {
        let x = self.tau(t);
        [0, 1, 2].map(|i| cheb_eval(&self.coef[i], x))
    }

    /// Velocity from the same polynomial; the chain rule factor is 2 / span.
    pub fn velocity(&self, t: f64) -> [f64; 3] {
        let x = self.tau(t);
        let s = 2.0 / (self.t1 - self.t0);
        [0, 1, 2].map(|i| s * cheb_eval(&cheb_deriv(&self.coef[i]), x))
    }

    /// Max position (km) and velocity (km/s) error on a dense interior grid
    /// that deliberately avoids the fit nodes.
    pub fn max_error(&self, orb: &Orbit, samples: usize) -> (f64, f64) {
        let (mut ep, mut ev) = (0.0f64, 0.0f64);
        for j in 0..=samples {
            let t = self.t0 + (self.t1 - self.t0) * j as f64 / samples as f64;
            let (r, v) = orb.state(t);
            ep = ep.max(dist(&self.position(t), &r));
            ev = ev.max(dist(&self.velocity(t), &v));
        }
        (ep, ev)
    }
}

pub fn dist(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

pub const CHECK_SAMPLES: usize = 257;

/// Equal-length segments over one period: smallest count meeting `tol` (km).
pub fn uniform_layout(orb: &Orbit, n: usize, tol: f64, kmax: usize) -> Option<Vec<Segment>> {
    let p = orb.period();
    (1..=kmax).find_map(|k| {
        let segs: Vec<Segment> = (0..k)
            .map(|i| Segment::fit(orb, p * i as f64 / k as f64, p * (i + 1) as f64 / k as f64, n))
            .collect();
        segs.iter()
            .all(|s| s.max_error(orb, CHECK_SAMPLES).0 <= tol)
            .then_some(segs)
    })
}

/// Recursive bisection: split any segment whose checked error exceeds `tol`.
/// Starts from two half-period segments split at perigee/apogee.
pub fn adaptive_layout(orb: &Orbit, n: usize, tol: f64) -> Vec<Segment> {
    let p = orb.period();
    let mut out = Vec::new();
    let mut stack = vec![(0.5 * p, p), (0.0, 0.5 * p)];
    while let Some((t0, t1)) = stack.pop() {
        let s = Segment::fit(orb, t0, t1, n);
        if s.max_error(orb, CHECK_SAMPLES).0 <= tol || t1 - t0 < 1e-3 {
            out.push(s);
        } else {
            let tm = 0.5 * (t0 + t1);
            stack.push((tm, t1));
            stack.push((t0, tm));
        }
    }
    out
}

/// A-priori layout from the singularity geometry alone: march from perigee
/// and end each segment where its Bernstein parameter drops to `rho_star`,
/// so every segment converges at the same rate. Short near perigee, long
/// near apogee, and no power-of-two restriction.
pub fn equal_rho_layout(orb: &Orbit, n: usize, rho_star: f64) -> Vec<Segment> {
    let p = orb.period();
    let mut out = Vec::new();
    let mut t0 = 0.0;
    while t0 < p - 1e-9 {
        let mut t1 = p;
        if segment_rho(orb, t0, p) < rho_star {
            // rho falls monotonically as the segment grows: bisect on the end.
            let (mut lo, mut hi) = (t0, p);
            for _ in 0..60 {
                let mid = 0.5 * (lo + hi);
                if segment_rho(orb, t0, mid) >= rho_star { lo = mid } else { hi = mid }
            }
            t1 = lo;
        }
        out.push(Segment::fit(orb, t0, t1, n));
        t0 = t1;
    }
    out
}

/// Smallest rho* (fewest segments) whose equal-rho layout meets `tol`.
/// Only this outer acceptance test touches the truth model.
pub fn tuned_equal_rho_layout(orb: &Orbit, n: usize, tol: f64) -> (f64, Vec<Segment>) {
    let ok = |r: f64| {
        let segs = equal_rho_layout(orb, n, r);
        let pass = segs.iter().all(|s| s.max_error(orb, CHECK_SAMPLES).0 <= tol);
        (pass, segs)
    };
    let (mut lo, mut hi) = (1.05f64, 200.0f64);
    assert!(ok(hi).0, "tolerance unreachable at this degree");
    for _ in 0..30 {
        let mid = (lo * hi).sqrt();
        if ok(mid).0 { hi = mid } else { lo = mid }
    }
    (hi, ok(hi).1)
}

// --------------------------------------------------- Hermite baseline

/// Cubic Hermite through position and velocity at uniform knots: the usual
/// "store state vectors every dt" ephemeris. Returns max position error.
pub fn hermite_error(orb: &Orbit, knots: usize, samples_per_span: usize) -> f64 {
    let p = orb.period();
    let h = p / knots as f64;
    let mut emax = 0.0f64;
    for i in 0..knots {
        let (t0, t1) = (i as f64 * h, (i + 1) as f64 * h);
        let ((r0, v0), (r1, v1)) = (orb.state(t0), orb.state(t1));
        for j in 1..samples_per_span {
            let s = j as f64 / samples_per_span as f64;
            let (h00, h10) = (2.0 * s.powi(3) - 3.0 * s * s + 1.0, s.powi(3) - 2.0 * s * s + s);
            let (h01, h11) = (-2.0 * s.powi(3) + 3.0 * s * s, s.powi(3) - s * s);
            let r: [f64; 3] = [0, 1, 2].map(|a| h00 * r0[a] + h10 * h * v0[a] + h01 * r1[a] + h11 * h * v1[a]);
            emax = emax.max(dist(&r, &orb.state(t0 + s * h).0));
        }
    }
    emax
}

pub fn hermite_knots_for(orb: &Orbit, tol: f64) -> usize {
    // Error falls as knots^-4; bracket by doubling then bisect on the count.
    let mut hi = 8;
    while hermite_error(orb, hi, 32) > tol {
        hi *= 2;
    }
    let mut lo = hi / 2;
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if hermite_error(orb, mid, 32) > tol {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    hi
}

/// Upper-envelope decay rate of |c_k|: least-squares slope of ln|c_k| over
/// the running maxima (from the tail) in [k_lo, k_hi], with the k^-3/2
/// algebraic factor of a square-root branch point removed. Returns the
/// estimated rho = exp(-slope).
pub fn measured_rho(c: &[f64], k_lo: usize, k_hi: usize) -> f64 {
    let mut pts = Vec::new();
    for k in k_lo..=k_hi {
        let env = c[k..=(k + 3).min(c.len() - 1)].iter().fold(0.0f64, |m, v| m.max(v.abs()));
        pts.push((k as f64, env.ln() + 1.5 * (k as f64).ln()));
    }
    let nn = pts.len() as f64;
    let (sx, sy) = pts.iter().fold((0.0, 0.0), |(a, b), p| (a + p.0, b + p.1));
    let (mx, my) = (sx / nn, sy / nn);
    let (num, den) = pts
        .iter()
        .fold((0.0, 0.0), |(a, b), p| (a + (p.0 - mx) * (p.1 - my), b + (p.0 - mx).powi(2)));
    (-num / den).exp()
}
