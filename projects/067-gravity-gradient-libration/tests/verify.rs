// Astron Built this AP
// Project 067 | Attitude control | Gravity-gradient stability regions and eccentric-orbit pitch libration
// Author: Aditya (aditya-astron)
// Date: 2026-11-27 IST
//
// The nonlinear 3D model is the reference. Each reduced model (linear
// roll/yaw, linear pitch, eccentric planar pitch, Floquet theory) is checked
// against it or against a closed form before the design study relies on it.

use gg_libration::linear::{self, Region};
use gg_libration::ode::{char_poly, poly_from_roots, rk4, transition_matrix};
use gg_libration::pitch;
use gg_libration::rigid::{dcm, mv, quat_from_euler, Body};
use std::f64::consts::PI;

fn rot(axis: usize, a: f64) -> [[f64; 3]; 3] {
    // passive (frame) rotation
    let (c, s) = (a.cos(), a.sin());
    match axis {
        1 => [[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]],
        2 => [[c, 0.0, -s], [0.0, 1.0, 0.0], [s, 0.0, c]],
        _ => [[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]],
    }
}

fn mm(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

fn propagate(b: &Body, y: &mut Vec<f64>, t1: f64, steps: usize) {
    let f = |t: f64, s: &[f64], d: &mut [f64]| b.rhs(t, s, d);
    rk4(&f, 0.0, t1, y, steps, |_, _| {});
}

#[test]
fn quaternion_matches_321_dcm() {
    let (r, p, y) = (0.3, -0.7, 1.1);
    let c = dcm(&quat_from_euler(r, p, y));
    let reference = mm(&mm(&rot(1, r), &rot(2, p)), &rot(3, y));
    for i in 0..3 {
        for j in 0..3 {
            assert!((c[i][j] - reference[i][j]).abs() < 1e-14);
        }
    }
}

#[test]
fn quaternion_rate_gives_dcm_rate() {
    // C' = -[w_bo x] C, checked by central difference of the integrated state
    let b = Body::from_k(0.6, 0.3, 0.0);
    let y0 = b.state(0.0, (0.2, -0.4, 0.9), (0.05, -0.03, 0.07));
    let mut dy = vec![0.0; 8];
    b.rhs(0.0, &y0, &mut dy);
    let h = 1e-6;
    let qp: Vec<f64> = (0..4).map(|k| y0[k] + h * dy[k]).collect();
    let qm: Vec<f64> = (0..4).map(|k| y0[k] - h * dy[k]).collect();
    let (cp, cm, c) = (dcm(&qp), dcm(&qm), dcm(&y0[0..4]));
    let w_oi = mv(&c, &[0.0, -1.0, 0.0]);
    let w = [y0[4] - w_oi[0], y0[5] - w_oi[1], y0[6] - w_oi[2]];
    for j in 0..3 {
        let col = [c[0][j], c[1][j], c[2][j]];
        let wxc = [
            w[1] * col[2] - w[2] * col[1],
            w[2] * col[0] - w[0] * col[2],
            w[0] * col[1] - w[1] * col[0],
        ];
        for i in 0..3 {
            let fd = (cp[i][j] - cm[i][j]) / (2.0 * h);
            assert!((fd + wxc[i]).abs() < 1e-8, "C' mismatch at {},{}", i, j);
        }
    }
}

#[test]
fn jacobi_integral_conserved_circular() {
    // 10 deg tumble-free libration, 10 orbits; RK4 error must fall ~16x per halving
    let b = Body::from_k(0.7, 0.35, 0.0);
    let y0 = b.state(0.0, (0.17, 0.17, -0.17), (0.02, 0.0, -0.01));
    let h0 = b.jacobi(&y0);
    let mut errs = vec![];
    for &n in &[4000usize, 8000] {
        let mut y = y0.clone();
        propagate(&b, &mut y, 20.0 * PI, n);
        errs.push((b.jacobi(&y) - h0).abs());
    }
    assert!(errs[1] < 1e-9, "Jacobi drift {:e}", errs[1]);
    assert!(errs[0] / errs[1] > 10.0, "not fourth order: {:?}", errs);
}

#[test]
fn linear_roll_yaw_matches_nonlinear() {
    // The potential is even in the angles, so there is no quadratic term:
    // the linear model's relative error must fall as a^2.
    for &(k1, k3) in &[(0.8, 0.4), (0.3, 0.1), (-0.05, -0.2)] {
        let b = Body::from_k(k1, k3, 0.0);
        let mut rel = vec![];
        for &a in &[1e-3, 5e-4] {
            let mut y = b.state(0.0, (a, 0.0, -0.7 * a), (0.0, 0.0, 0.3 * a));
            propagate(&b, &mut y, 6.0 * PI, 9000);
            let (r, _p, yw) = b.angles(&y);
            let mut x = vec![a, -0.7 * a, 0.0, 0.3 * a];
            rk4(
                &linear::roll_yaw_rhs(k1, k3),
                0.0,
                6.0 * PI,
                &mut x,
                9000,
                |_, _| {},
            );
            rel.push(((r - x[0]).abs() + (yw - x[1]).abs()) / a);
        }
        assert!(rel[0] < 1e-4, "k=({},{}) rel {:e}", k1, k3, rel[0]);
        let ratio = rel[0] / rel[1];
        assert!(
            ratio > 3.6 && ratio < 4.4,
            "relative error not O(a^2): {:?}",
            rel
        );
    }
}

#[test]
fn linear_pitch_matches_nonlinear() {
    let b = Body::from_k(0.8, 0.4, 0.0);
    let k2 = b.k2();
    let a = 1e-4;
    let mut y = b.state(0.0, (0.0, a, 0.0), (0.0, 0.0, 0.0));
    propagate(&b, &mut y, 3.0 * PI, 6000);
    let w = linear::pitch_freq(k2).unwrap();
    let expect = a * (w * 3.0 * PI).cos();
    assert!((b.angles(&y).1 - expect).abs() < 1e-10);
}

#[test]
fn characteristic_exponents_vs_transition_matrix() {
    // eigenvalues of Phi(T) must be exp(s_i T) for the analytic s_i
    for &(k1, k3) in &[
        (0.8, 0.4),
        (0.9, 0.85),
        (-0.05, -0.2),
        (0.4, -0.3),
        (-0.5, -0.5),
    ] {
        let t = 1.3;
        let phi = transition_matrix(&linear::roll_yaw_rhs(k1, k3), 4, 0.0, t, 4000);
        let cp = char_poly(&phi, 4);
        let ex: Vec<(f64, f64)> = linear::roll_yaw_roots(k1, k3)
            .iter()
            .map(|&(re, im)| {
                let m = (re * t).exp();
                (m * (im * t).cos(), m * (im * t).sin())
            })
            .collect();
        let pr = poly_from_roots(&ex);
        for k in 0..5 {
            let scale = 1.0 + cp[k].abs();
            assert!(
                (cp[k] - pr[k].0).abs() < 1e-9 * scale,
                "k=({},{}) coeff {}",
                k1,
                k3,
                k
            );
            assert!(pr[k].1.abs() < 1e-9 * scale);
        }
    }
}

#[test]
fn reference_frequencies() {
    // thin boom (I3 -> 0, k1 = 1): roll/yaw at 2 w0 and sqrt(k3) w0,
    // pitch at sqrt(3) w0 when also k2 = 1
    for &k3 in &[0.1, 0.5, 0.9] {
        let (hi, lo) = linear::roll_yaw_freqs(1.0, k3).unwrap();
        assert!((hi - 2.0).abs() < 1e-12 && (lo - k3.sqrt()).abs() < 1e-12);
    }
    assert!((linear::pitch_freq(1.0).unwrap() - 3f64.sqrt()).abs() < 1e-15);
    assert!((linear::k2_from(0.8, 0.4) - Body::from_k(0.8, 0.4, 0.0).k2()).abs() < 1e-14);
}

#[test]
fn faddeev_leverrier_known_matrix() {
    // [[2,1],[1,3]]: lambda^2 - 5 lambda + 5
    let c = char_poly(&[2.0, 1.0, 1.0, 3.0], 2);
    assert!((c[1] + 5.0).abs() < 1e-14 && (c[2] - 5.0).abs() < 1e-14);
}

#[test]
fn region_classification() {
    assert_eq!(linear::classify(0.8, 0.4), Region::Lagrange);
    assert_eq!(linear::classify(0.4, 0.8), Region::Unstable); // pitch: I1 < I3
    assert_eq!(linear::classify(-0.05, -0.2), Region::DeBraDelp);
    assert_eq!(linear::classify(-0.5, -0.5), Region::Unstable); // b^2 < 4c
    assert_eq!(linear::classify(0.4, -0.3), Region::Unstable); // c < 0
                                                               // Lagrange area is the triangle k1 > k3 > 0, 1/8 of the k-square
    let m = 801;
    let mut n = 0;
    for a in 0..m {
        for b in 0..m {
            let k1 = -1.0 + 2.0 * (a as f64 + 0.5) / m as f64;
            let k3 = -1.0 + 2.0 * (b as f64 + 0.5) / m as f64;
            if linear::classify(k1, k3) == Region::Lagrange {
                n += 1;
                assert!(k1 > k3 && k3 > 0.0);
            }
        }
    }
    assert!((n as f64 / (m * m) as f64 - 0.125).abs() < 2e-3);
}

#[test]
fn debra_delp_is_gyroscopic_not_potential() {
    // DeBra-Delp: bounded in the full model, but the Jacobi potential is a
    // saddle (a pure roll offset lowers it); Lagrange is a true minimum.
    let dd = Body::from_k(-0.05, -0.2, 0.0);
    let lg = Body::from_k(0.8, 0.4, 0.0);
    let still =
        |b: &Body, r: f64, p: f64, y: f64| b.jacobi(&b.state(0.0, (r, p, y), (0.0, 0.0, 0.0)));
    let d = 1e-3;
    assert!(still(&dd, d, 0.0, 0.0) < still(&dd, 0.0, 0.0, 0.0));
    for &(r, p, y) in &[(d, 0.0, 0.0), (0.0, d, 0.0), (0.0, 0.0, d)] {
        assert!(still(&lg, r, p, y) > still(&lg, 0.0, 0.0, 0.0));
    }
    let a = 1e-3;
    let mut y = dd.state(0.0, (a, 0.0, a), (0.0, 0.0, 0.0));
    let mut peak = 0.0f64;
    let f = |t: f64, s: &[f64], ds: &mut [f64]| dd.rhs(t, s, ds);
    rk4(&f, 0.0, 40.0 * PI, &mut y, 40000, |_, s| {
        let (r, p, yw) = dd.angles(s);
        peak = peak.max(r.abs().max(p.abs()).max(yw.abs()));
    });
    assert!(peak < 10.0 * a, "DeBra-Delp grew to {:e}", peak);

    // and an unstable point does grow
    let un = Body::from_k(0.4, -0.3, 0.0);
    let mut y = un.state(0.0, (a, 0.0, a), (0.0, 0.0, 0.0));
    propagate(&un, &mut y, 10.0 * PI, 10000);
    let (r, _, yw) = un.angles(&y);
    assert!(r.abs().max(yw.abs()) > 100.0 * a);
}

#[test]
fn planar_eccentric_equation_matches_3d() {
    for &(e, th0, thd) in &[(0.1, 0.05, 0.02), (0.3, -0.4, 0.1)] {
        let b = Body::from_k(0.9, 0.0, e);
        let mut y = b.state(0.0, (0.0, th0, 0.0), (0.0, thd * b.nu_dot(0.0), 0.0));
        // one anomalistic period is 2pi in time, and nu goes 0 -> 2pi
        propagate(&b, &mut y, 2.0 * PI, 40000);
        let mut x = vec![th0, thd];
        rk4(
            &pitch::rhs_nonlinear(b.k2(), e),
            0.0,
            2.0 * PI,
            &mut x,
            40000,
            |_, _| {},
        );
        assert!((y[7] - 2.0 * PI).abs() < 1e-10);
        assert!(
            (b.angles(&y).1 - x[0]).abs() < 1e-9,
            "e={} {} vs {}",
            e,
            b.angles(&y).1,
            x[0]
        );
    }
}

#[test]
fn monodromy_properties() {
    // det = 1 by Liouville (zero mean damping)
    for &(k2, e) in &[(0.2, 0.1), (0.9, 0.3), (0.08, 0.05)] {
        let m = pitch::monodromy(k2, e, 4000);
        assert!((m[0] * m[3] - m[1] * m[2] - 1.0).abs() < 1e-11);
    }
    // circular: tr = 2 cos(2 pi sqrt(3 k2))
    for &k2 in &[0.1, 0.5, 0.9] {
        let tr = pitch::trace(k2, 0.0);
        assert!((tr - 2.0 * (2.0 * PI * (3.0f64 * k2).sqrt()).cos()).abs() < 1e-11);
    }
    // k2 = 1/3: coexistence, M = I at any eccentricity
    for &e in &[0.01, 0.1, 0.4] {
        let m = pitch::monodromy(1.0 / 3.0, e, 4000);
        assert!((m[0] - 1.0).abs() + m[1].abs() + m[2].abs() + (m[3] - 1.0).abs() < 1e-11);
    }
}

#[test]
fn first_tongue_matches_mathieu_prediction() {
    // edges 1/12 -+ e/8 with an O(e^2) error
    let mut prev = None;
    for &e in &[0.04, 0.02, 0.01] {
        let (lo, hi) = pitch::tongue(1, e).unwrap();
        let err = (lo - (1.0 / 12.0 - e / 8.0)).abs() + (hi - (1.0 / 12.0 + e / 8.0)).abs();
        assert!(err < 0.2 * e * e + 1e-9, "e={} err {:e}", e, err);
        assert!(((hi - lo) / (e / 4.0) - 1.0).abs() < 0.05);
        if let Some(p) = prev {
            let r: f64 = p / err;
            assert!(r > 3.0 && r < 5.0, "error ratio {}", r);
        }
        prev = Some(err);
    }
}

#[test]
fn third_tongue_width_scales_as_e_cubed() {
    let w = |e| pitch::tongue(3, e).map(|(a, b)| b - a).unwrap();
    let r = w(0.2) / w(0.1);
    assert!((r / 8.0 - 1.0).abs() < 0.05, "ratio {}", r);
    assert!(pitch::tongue(2, 0.2).is_none());
}

#[test]
fn growth_inside_tongue_matches_multiplier() {
    let e = 0.05;
    let (lo, hi) = pitch::tongue(1, e).unwrap();
    let k2 = 0.5 * (lo + hi);
    let lam = pitch::max_multiplier(pitch::trace(k2, e));
    assert!(lam > 1.05);
    // x_n = a lambda^n v1 + b lambda^-n v2, so after enough orbits the
    // orbit-to-orbit growth is lambda
    let mut x = vec![1e-6, 0.0];
    let mut norms = vec![];
    for _ in 0..60 {
        rk4(
            &pitch::rhs_linear(k2, e),
            0.0,
            2.0 * PI,
            &mut x,
            4000,
            |_, _| {},
        );
        norms.push((x[0] * x[0] + x[1] * x[1]).sqrt());
    }
    let rate = (norms[59] / norms[58]).ln();
    assert!(
        (rate / lam.ln() - 1.0).abs() < 1e-6,
        "{} vs {}",
        rate,
        lam.ln()
    );
}

#[test]
fn periodic_libration_small_e_matches_linear() {
    let k2 = 0.9;
    let mut prev = None;
    for &e in &[0.02, 0.01, 0.005] {
        let a = pitch::forced_amplitude_linear(k2, e);
        let (x0, m) = pitch::periodic_orbit(k2, e, [0.0, a]).unwrap();
        let pk = pitch::peak_over_orbit(k2, e, x0);
        // the peak of the exact periodic libration differs from |A| at O(e^2)
        let rel = (pk / a.abs() - 1.0).abs();
        assert!(rel < 5.0 * e * e, "e={} rel {}", e, rel);
        assert!((m[0] + m[3]).abs() < 2.0);
        if let Some(p) = prev {
            let r: f64 = p / rel;
            assert!(r > 3.6 && r < 4.4, "not O(e^2): {}", r);
        }
        prev = Some(rel);
    }
}

#[test]
fn periodic_orbit_closes_in_3d() {
    let e = 0.05;
    let b = Body::from_k(0.9, 0.2, e);
    let k2 = b.k2();
    let (x0, _) =
        pitch::periodic_orbit(k2, e, [0.0, pitch::forced_amplitude_linear(k2, e)]).unwrap();
    let mut y = b.state(0.0, (0.0, x0[0], 0.0), (0.0, x0[1] * b.nu_dot(0.0), 0.0));
    propagate(&b, &mut y, 10.0 * PI, 30000);
    assert!((b.angles(&y).1 - x0[0]).abs() < 1e-9);
}
