// Astron Built this AP
// Project 029 | Attitude control | B-dot magnetic detumbling of a 3U CubeSat
// Author: Aditya (aditya-astron)
// Date: 2026-11-11 IST

use bdot_detumble::*;

const D2R: f64 = std::f64::consts::PI / 180.0;

fn close(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol
}

#[test]
fn dipole_matches_closed_form() {
    let d = Dipole::aligned(3.0e-5);
    // Equator at the surface: B = B0 pointing north (+z).
    let b = d.field([R_EARTH, 0.0, 0.0], 0.0);
    assert!(close(b[2], 3.0e-5, 1e-18) && b[0].abs() < 1e-18);
    // North pole: |B| = 2 B0, pointing down into the Earth.
    let b = d.field([0.0, 0.0, R_EARTH], 0.0);
    assert!(close(b[2], -6.0e-5, 1e-18));
    // Inverse cube with altitude.
    let b = d.field([2.0 * R_EARTH, 0.0, 0.0], 0.0);
    assert!(close(b[2], 3.0e-5 / 8.0, 1e-18));

    // IGRF-2020 dipole: B0 ~ 29.8 uT, axis tilted ~9.4 deg (pole near 80.7 N).
    let g = Dipole::igrf2020();
    assert!(close(g.b0 * 1e9, 29806.0, 1.0));
    let tilt = (-g.axis_ecef[2]).acos() / D2R;
    assert!(close(tilt, 9.4, 0.1), "tilt {tilt}");
}

#[test]
fn dipole_is_divergence_and_curl_free() {
    let d = Dipole::igrf2020();
    let r0 = [4.1e6, -3.3e6, 4.4e6];
    let h = 10.0;
    let mut grad = [[0.0; 3]; 3]; // grad[i][j] = dB_i / dx_j
    for j in 0..3 {
        let mut rp = r0;
        let mut rm = r0;
        rp[j] += h;
        rm[j] -= h;
        let (bp, bm) = (d.field(rp, 1234.0), d.field(rm, 1234.0));
        for i in 0..3 {
            grad[i][j] = (bp[i] - bm[i]) / (2.0 * h);
        }
    }
    let scale_g = 3.0 * norm(d.field(r0, 1234.0)) / norm(r0);
    let div = grad[0][0] + grad[1][1] + grad[2][2];
    assert!(div.abs() < 1e-7 * scale_g);
    for (i, j) in [(0, 1), (0, 2), (1, 2)] {
        assert!((grad[i][j] - grad[j][i]).abs() < 1e-7 * scale_g);
    }
}

#[test]
fn circular_orbit_period_and_geometry() {
    let o = Orbit { radius: R_EARTH + 500e3, incl: 97.4 * D2R, raan: 0.3, u0: 0.2 };
    // Hand value: 2 pi sqrt(6871.2e3^3 / 3.986004418e14) = 2 pi x 902.15 = 5668.4 s
    assert!(close(o.period(), 5668.4, 0.1), "{}", o.period());
    for t in [0.0, 1000.0, 3000.0] {
        assert!(close(norm(o.position(t)), o.radius, 1e-6));
    }
    let p = o.position(0.5 * o.period());
    assert!(close(norm(add(p, o.position(0.0))), 0.0, 1e-6));
}

#[test]
fn jacobi_eigenvalues() {
    let ev = sym_eigenvalues([[2.0, 1.0, 0.0], [1.0, 2.0, 1.0], [0.0, 1.0, 2.0]]);
    let s2 = 2f64.sqrt();
    assert!(close(ev[0], 2.0 - s2, 1e-13) && close(ev[1], 2.0, 1e-13) && close(ev[2], 2.0 + s2, 1e-13));
    let ev = sym_eigenvalues([[4.0, -2.0, 1.0], [-2.0, 3.0, 0.5], [1.0, 0.5, 6.0]]);
    // trace and determinant invariants
    assert!(close(ev.iter().sum::<f64>(), 13.0, 1e-12));
    let det = 4.0 * (3.0 * 6.0 - 0.25) + 2.0 * (-2.0 * 6.0 - 0.5) + 1.0 * (-1.0 - 3.0);
    assert!(close(ev[0] * ev[1] * ev[2], det, 1e-11));
}

#[test]
fn saturation_keeps_direction() {
    let m = saturate([0.5, -0.1, 0.3], 0.2);
    assert!(close(m[0], 0.2, 1e-15) && close(m[1], -0.04, 1e-15) && close(m[2], 0.12, 1e-15));
    assert_eq!(saturate([0.1, 0.0, -0.15], 0.2), [0.1, 0.0, -0.15]);
}

#[test]
fn averaged_field_direction_closed_forms() {
    // Equatorial orbit: b is always along the axis, spin about it is undamped.
    let ev = sym_eigenvalues(averaged_bbt(0.0, 360));
    assert!(close(ev[2], 1.0, 1e-14) && ev[0].abs() < 1e-14 && ev[1].abs() < 1e-14);

    // Polar orbit (orbit plane x-z): |B|^2 ~ 1 + 3 sin^2 u, and with
    // <1/(1+3 sin^2 u)> = 1/2 the averages come out as <bx^2> = <bz^2> = 1/2,
    // <by^2> = 0, off-diagonals 0 by symmetry.
    let m = averaged_bbt(90.0 * D2R, 360);
    assert!(close(m[0][0], 0.5, 1e-12) && close(m[2][2], 0.5, 1e-12) && m[1][1].abs() < 1e-12);
    assert!(m[0][2].abs() < 1e-12 && m[0][1].abs() < 1e-12);

    for xi in [10.0, 35.0, 60.0] {
        let m = averaged_bbt(xi * D2R, 360);
        assert!(close(m[0][0] + m[1][1] + m[2][2], 1.0, 1e-13));
        // geometric convergence of the periodic midpoint rule
        let m2 = averaged_bbt(xi * D2R, 720);
        assert!((m[1][1] - m2[1][1]).abs() < 1e-13);
    }
}

#[test]
fn torque_free_motion_conserves_energy_and_momentum() {
    let s = Scenario {
        sc: Spacecraft { inertia: [0.0405, 0.031, 0.0067], m_max: 0.2 },
        orbit: Orbit { radius: R_EARTH + 500e3, incl: 1.0, raan: 0.0, u0: 0.0 },
        field: Dipole::igrf2020(),
        law: Law::Ideal { k: 0.0 },
        w0: [0.05, 0.12, -0.2],
        q0: [1.0, 0.0, 0.0, 0.0],
        dt: 0.02,
        seed: 1,
    };
    let log = simulate(&s, 600.0, 600.0);
    let (a, b) = (log.first().unwrap(), log.last().unwrap());
    assert!((b.energy / a.energy - 1.0).abs() < 1e-9);
    // |H| is frame independent
    let j = s.sc.inertia;
    let h = |w: V3| norm([j[0] * w[0], j[1] * w[1], j[2] * w[2]]);
    assert!((h(b.w) / h(a.w) - 1.0).abs() < 1e-9);
}

#[test]
fn ideal_law_is_dissipative_every_step() {
    // dE/dt = w . tau = -k |w x b|^2 <= 0
    let s = Scenario {
        sc: Spacecraft { inertia: [0.0405, 0.0405, 0.0067], m_max: 0.2 },
        orbit: Orbit { radius: R_EARTH + 500e3, incl: 51.6 * D2R, raan: 0.0, u0: 0.0 },
        field: Dipole::igrf2020(),
        law: Law::Ideal { k: 3e-5 },
        w0: [0.1, -0.15, 0.4],
        q0: [0.5, 0.5, 0.5, 0.5],
        dt: 0.05,
        seed: 1,
    };
    let log = simulate(&s, 3000.0, 0.05);
    for p in log.windows(2) {
        assert!(p[1].energy <= p[0].energy * (1.0 + 1e-12));
    }
    assert!(log.last().unwrap().energy < 0.2 * log[0].energy);
}

/// Floquet-style check of the averaging analysis: isotropic body, polar orbit
/// in an aligned dipole, small gain. Theory says the in-plane spin decays at
/// 0.5 k/J; measure the log slope of |w| between whole orbits.
fn measured_rate(eps: f64) -> (f64, f64) {
    let j = 0.02;
    let orbit = Orbit { radius: R_EARTH + 500e3, incl: 90.0 * D2R, raan: 0.0, u0: 0.0 };
    let n = orbit.mean_motion();
    let k = eps * n * j;
    let s = Scenario {
        sc: Spacecraft { inertia: [j, j, j], m_max: 1.0 },
        orbit,
        field: Dipole::aligned(3.0e-5),
        law: Law::Ideal { k },
        w0: [0.01, 0.003, 0.008],
        q0: [0.8, 0.0, 0.6, 0.0],
        dt: 0.5,
        seed: 1,
    };
    let tp = orbit.period();
    let (o1, o2) = (12.0, 24.0);
    let log = simulate(&s, o2 * tp, tp / 4.0);
    let at = |o: f64| norm(log[(o * 4.0).round() as usize].w);
    let measured = (at(o1) / at(o2)).ln() / ((o2 - o1) * tp);
    (measured, averaged_slowest_rate(k, j, 90.0 * D2R))
}

#[test]
fn averaging_theory_predicts_decay_rate() {
    let (m1, p1) = measured_rate(0.10);
    let (m2, p2) = measured_rate(0.05);
    let e1 = (m1 / p1 - 1.0).abs();
    let e2 = (m2 / p2 - 1.0).abs();
    println!("eps 0.10: measured {m1:.4e}, averaged {p1:.4e}, err {:.2} %", 100.0 * e1);
    println!("eps 0.05: measured {m2:.4e}, averaged {p2:.4e}, err {:.2} %", 100.0 * e2);
    assert!(e1 < 0.03 && e2 < 0.03);
}

#[test]
fn sampled_law_tracks_ideal_law_at_fast_spin() {
    let base = Scenario {
        sc: Spacecraft { inertia: [0.0405, 0.0405, 0.0067], m_max: 100.0 },
        orbit: Orbit { radius: R_EARTH + 500e3, incl: 97.4 * D2R, raan: 0.5, u0: 0.0 },
        field: Dipole::igrf2020(),
        law: Law::Ideal { k: 3e-5 },
        w0: [10.0 * D2R, -8.0 * D2R, 6.0 * D2R],
        q0: [1.0, 0.0, 0.0, 0.0],
        dt: 0.01,
        seed: 1,
    };
    let mut sampled = base;
    sampled.law = Law::Sampled { k: 3e-5, period: 0.1, quiet: 0.0, noise: 0.0 };
    let t = 1000.0;
    let wi = norm(simulate(&base, t, t).last().unwrap().w);
    let ws = norm(simulate(&sampled, t, t).last().unwrap().w);
    println!("|w| after {t} s: ideal {:.4} deg/s, sampled {:.4} deg/s", wi / D2R, ws / D2R);
    assert!((ws / wi - 1.0).abs() < 0.03);
    assert!(wi < 0.6 * norm(base.w0));
}

#[test]
fn reference_detumble_reaches_twice_orbit_rate() {
    let orbit = Orbit { radius: R_EARTH + 500e3, incl: 97.4 * D2R, raan: 30.0 * D2R, u0: 0.0 };
    let sc = Spacecraft { inertia: [0.0405, 0.0405, 0.0067], m_max: 0.2 };
    let n = orbit.mean_motion();
    let k = avanzini_gain(n, (180.0 - 97.4) * D2R, sc.j_min());
    let s = Scenario {
        sc,
        orbit,
        field: Dipole::igrf2020(),
        law: Law::Sampled { k, period: 1.0, quiet: 0.05, noise: 50e-9 },
        w0: [10.0 * D2R, -8.0 * D2R, 6.0 * D2R],
        q0: [0.9, 0.1, -0.3, 0.3],
        dt: 0.05,
        seed: 2029,
    };
    let tp = orbit.period();
    let log = simulate(&s, 6.0 * tp, 10.0);
    let t_settle = settle_time(&log, 0.5 * D2R).expect("did not detumble");
    assert!(t_settle < 3.0 * tp);
    assert!(log.iter().all(|p| p.dipole <= sc.m_max * (1.0 + 1e-12)));
    // B-dot leaves the body turning with the field, which in a near-polar
    // orbit sweeps round about twice per orbit.
    let tail: Vec<f64> = log.iter().filter(|p| p.t > 5.0 * tp).map(|p| norm(p.w)).collect();
    let mean = tail.iter().sum::<f64>() / tail.len() as f64;
    assert!(mean > 1.5 * n && mean < 2.5 * n, "residual {} n", mean / n);
}

#[test]
fn settle_time_semantics() {
    let mk = |t: f64, w: f64| Sample { t, w: [w, 0.0, 0.0], energy: 0.0, dipole: 0.0 };
    let log = vec![mk(0.0, 5.0), mk(1.0, 0.5), mk(2.0, 2.0), mk(3.0, 0.4), mk(4.0, 0.3)];
    assert_eq!(settle_time(&log, 1.0), Some(3.0));
    assert_eq!(settle_time(&log, 0.1), None);
}
