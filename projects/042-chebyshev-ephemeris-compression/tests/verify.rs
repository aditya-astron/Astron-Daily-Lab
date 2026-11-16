// Astron Built this AP
// Project 042 | Numerical methods / Navigation | Chebyshev ephemeris compression for an eccentric orbit
// Author: Aditya (aditya-astron)
// Date: 2026-11-16 IST

use cheb_ephemeris::complex::C64;
use cheb_ephemeris::*;
use std::f64::consts::PI;

fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol
}

fn molniya() -> Orbit {
    Orbit { a: 26_554.0, e: 0.74, inc: 63.4f64.to_radians(), raan: 0.7, argp: 1.5 * PI, mu: MU_EARTH }
}

// Jacobi-Anger / modified-Bessel generating functions:
//   exp(x)   = I0(1) + 2 sum I_k(1) T_k(x)
//   cos(a x) = J0(a) + 2 sum (-1)^k J_2k(a) T_2k(x)
// Bessel values from Abramowitz & Stegun tables 9.1 and 9.8 (also scipy.special).
#[test]
fn fit_reproduces_bessel_coefficients() {
    let c = cheb_fit(f64::exp, 20);
    let i_k = [1.2660658777520082, 0.5651591039924851, 0.1357476697670383, 0.0221684249243319, 0.0027371202210468];
    assert!(close(c[0], i_k[0], 1e-15));
    for k in 1..5 {
        assert!(close(c[k], 2.0 * i_k[k], 1e-15), "exp k={k}: {}", c[k]);
    }
    let c = cheb_fit(|x| (2.0 * x).cos(), 24);
    let (j0, j2, j4, j6) = (0.2238907791412357, 0.3528340286156377, 0.0339957198075684, 0.0012024289717899);
    assert!(close(c[0], j0, 1e-15));
    assert!(close(c[2], -2.0 * j2, 1e-15) && close(c[4], 2.0 * j4, 1e-15) && close(c[6], -2.0 * j6, 1e-15));
    assert!(c.iter().skip(1).step_by(2).all(|v| v.abs() < 1e-15), "odd terms of an even function");
}

#[test]
fn clenshaw_matches_trig_definition_and_derivative() {
    // Single T_k: Clenshaw against cos(k acos x); derivative against k U_{k-1}.
    for k in 0..12 {
        let mut c = vec![0.0; 12];
        c[k] = 1.0;
        let d = cheb_deriv(&c);
        for j in 0..=40 {
            let x = -0.999 + 1.998 * j as f64 / 40.0;
            let th = x.acos();
            assert!(close(cheb_eval(&c, x), (k as f64 * th).cos(), 1e-13));
            let u = if k == 0 { 0.0 } else { k as f64 * (k as f64 * th).sin() / th.sin() };
            assert!(close(cheb_eval(&d, x), u, 1e-11), "T_{k}' at {x}");
        }
    }
    // Interpolation is exact on polynomials of degree <= N.
    let p = |x: f64| 3.0 * x.powi(5) - x.powi(3) + 0.5 * x - 2.0;
    let c = cheb_fit(p, 5);
    for j in 0..=20 {
        let x = -1.0 + 0.1 * j as f64;
        assert!(close(cheb_eval(&c, x), p(x), 1e-13));
    }
}

#[test]
fn kepler_solver_and_two_body_invariants() {
    // Vallado, Fundamentals of Astrodynamics, Example 2-1: M = 235.4 deg, e = 0.4.
    let ea = solve_kepler(235.4f64.to_radians(), 0.4);
    assert!(close(ea.to_degrees(), 220.512_074_767_522, 1e-9));
    for &e in &[0.0, 0.3, 0.74, 0.95] {
        for j in 0..200 {
            let m = 2.0 * PI * j as f64 / 200.0;
            let ea = solve_kepler(m, e);
            assert!(close(ea - e * ea.sin(), m, 1e-12));
        }
    }
    let orb = molniya();
    let (r0, _) = orb.state(0.0);
    let rp = orb.a * (1.0 - orb.e);
    assert!(close(r0.iter().map(|x| x * x).sum::<f64>().sqrt(), rp, 1e-8));
    for j in 0..50 {
        let (r, v) = orb.state(orb.period() * j as f64 / 50.0);
        let rn = r.iter().map(|x| x * x).sum::<f64>().sqrt();
        let v2: f64 = v.iter().map(|x| x * x).sum();
        assert!(close(v2 / 2.0 - orb.mu / rn, -orb.mu / (2.0 * orb.a), 1e-10), "vis-viva");
    }
}

#[test]
fn bernstein_parameter_geometry() {
    // On [-1, 1] the ellipse degenerates: rho = 1.
    assert!(close(bernstein_rho(C64::new(0.3, 0.0)), 1.0, 1e-12));
    // Point i b on the minor axis: semi-minor (rho - 1/rho)/2 = b.
    for &b in &[0.1, 1.0, 5.0] {
        assert!(close(bernstein_rho(C64::new(0.0, b)), b + (1.0 + b * b).sqrt(), 1e-12));
    }
    // Point on the real axis beyond 1: semi-major (rho + 1/rho)/2 = x.
    assert!(close(bernstein_rho(C64::new(1.25, 0.0)), 2.0, 1e-12));
    // Branch point offset for e = 0.74 by hand: acosh(1/0.74) - sqrt(1 - 0.74^2).
    assert!(close(molniya().eta(), 0.142_882, 1e-6));
}

#[test]
fn coefficient_decay_matches_branch_point_prediction() {
    let orb = molniya();
    let n = orb.mean_motion();
    for &hm in &[0.15, 0.3, 0.6] {
        let (t0, t1) = (-hm / n, hm / n);
        let s = Segment::fit(&orb, t0, t1, 80);
        let pred = segment_rho(&orb, t0, t1);
        let c0 = s.coef[0][1].abs();
        let kmax = (1..s.coef[0].len()).filter(|&k| s.coef[0][k].abs() > 1e-11 * c0).max().unwrap();
        let meas = measured_rho(&s.coef[0], 6.max(kmax / 3), kmax - 3);
        assert!((meas / pred - 1.0).abs() < 0.01, "hm={hm}: pred {pred} meas {meas}");
    }
}

#[test]
fn layouts_meet_tolerance_on_an_independent_grid() {
    let orb = molniya();
    let tol = 1e-6;
    let (_, q) = tuned_equal_rho_layout(&orb, 13, tol);
    // Contiguous cover of exactly one period.
    assert!(close(q[0].t0, 0.0, 0.0) && close(q.last().unwrap().t1, orb.period(), 1e-6));
    assert!(q.windows(2).all(|w| w[0].t1 == w[1].t0));
    for s in &q {
        let (ep, ev) = s.max_error(&orb, 2003);
        assert!(ep < 1.05 * tol, "position {ep}");
        assert!(ev < 1e-6, "velocity {ev} km/s");
    }
    let u = uniform_layout(&orb, 13, tol, 4000).unwrap();
    let a = adaptive_layout(&orb, 13, tol);
    assert!(q.len() < a.len() && a.len() < u.len());
    assert!(q.len() * 3 * 14 * 10 < 6 * hermite_knots_for(&orb, tol), "Chebyshev should beat Hermite by >10x");
}

#[test]
fn hermite_error_is_fourth_order() {
    // Doubling the knots should cut the error by ~16 once in the asymptotic range.
    let orb = Orbit { e: 0.1, ..molniya() };
    let (e1, e2) = (hermite_error(&orb, 200, 32), hermite_error(&orb, 400, 32));
    let ratio = e1 / e2;
    assert!(ratio > 14.0 && ratio < 18.0, "ratio {ratio}");
}
