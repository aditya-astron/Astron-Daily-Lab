// Astron Built this AP
// Project 005 | Attitude determination | Wahba's problem: TRIAD, Davenport q-method and QUEST
// Author: Aditya (aditya-astron)
// Date: 2026-11-01 IST

use std::f64::consts::PI;
use wahba::*;

fn random_quat(rng: &mut Rng) -> Quat {
    quat_normalize([rng.gauss(), rng.gauss(), rng.gauss(), rng.gauss()])
}

fn noisy_set(rng: &mut Rng, a: &M3, n: usize, sigma: f64) -> Vec<Obs> {
    (0..n)
        .map(|_| {
            let r = rng.unit_vector();
            Obs { b: measure(rng, mat_vec(a, r), sigma), r, sigma }
        })
        .collect()
}

#[test]
fn dcm_matches_hand_rotation_about_z() {
    // Body frame rotated +90 deg about z from the reference frame:
    // reference x appears along body -y, reference y along body +x.
    let a = quat_to_dcm(quat_from_axis_angle([0.0, 0.0, 1.0], PI / 2.0));
    let expect = [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
    for i in 0..3 {
        for j in 0..3 {
            assert!((a[i][j] - expect[i][j]).abs() < 1e-15);
        }
    }
}

#[test]
fn quaternion_dcm_round_trip() {
    let mut rng = Rng::new(1);
    let mut qs: Vec<Quat> = (0..2000).map(|_| random_quat(&mut rng)).collect();
    // Shepperd branches: q4 = 0 and single dominant vector components.
    qs.push([1.0, 0.0, 0.0, 0.0]);
    qs.push([0.0, 1.0, 0.0, 0.0]);
    qs.push([0.0, 0.0, 1.0, 0.0]);
    qs.push(quat_from_axis_angle([1.0, 2.0, 3.0], PI));
    for q in qs {
        let a = quat_to_dcm(q);
        // Orthonormal with det +1.
        let i3 = mat_mul(&a, &transpose(&a));
        for r in 0..3 {
            for c in 0..3 {
                let e = if r == c { 1.0 } else { 0.0 };
                assert!((i3[r][c] - e).abs() < 1e-14);
            }
        }
        assert!((det3(&a) - 1.0).abs() < 1e-14);
        let q2 = dcm_to_quat(&a);
        let d = (0..4).map(|k| q[k] * q2[k]).sum::<f64>().abs();
        assert!((d - 1.0).abs() < 1e-13, "q and -q are the same attitude, nothing else is");
    }
}

#[test]
fn jacobi_recovers_known_spectrum() {
    // Q diag(4, -1, 2.5, 0.3) Q^T for an orthogonal Q built from a 4x4
    // Householder reflection.
    let v = [0.5, -0.5, 0.5, 0.5];
    let mut q = [[0.0; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            q[i][j] = if i == j { 1.0 } else { 0.0 } - 2.0 * v[i] * v[j];
        }
    }
    let lam = [4.0, -1.0, 2.5, 0.3];
    let mut m = [[0.0; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            m[i][j] = (0..4).map(|k| q[i][k] * lam[k] * q[j][k]).sum();
        }
    }
    let (vals, _) = jacobi_eig4(m);
    let mut got = vals.to_vec();
    got.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut want = lam.to_vec();
    want.sort_by(|a, b| a.partial_cmp(b).unwrap());
    for (g, w) in got.iter().zip(&want) {
        assert!((g - w).abs() < 1e-12);
    }
}

#[test]
fn noiseless_observations_are_recovered_exactly() {
    let mut rng = Rng::new(2);
    for _ in 0..500 {
        let a = quat_to_dcm(random_quat(&mut rng));
        let refs: Vec<V3> = (0..4).map(|_| rng.unit_vector()).collect();
        let obs: Vec<Obs> = refs.iter().map(|&r| Obs { b: mat_vec(&a, r), r, sigma: 1e-3 }).collect();
        let at = triad(obs[0].b, obs[1].b, obs[0].r, obs[1].r);
        assert!(error_angle(&at, &a) < 1e-12);
        assert!(error_angle(&quat_to_dcm(q_method(&obs).0), &a) < 1e-12);
        assert!(error_angle(&quat_to_dcm(quest(&obs).q), &a) < 1e-12);
    }
}

#[test]
fn quest_matches_q_method_and_wahba_identity() {
    let mut rng = Rng::new(3);
    for _ in 0..2000 {
        let a = quat_to_dcm(random_quat(&mut rng));
        let obs = noisy_set(&mut rng, &a, 5, 0.01);
        let qr = quest(&obs);
        let (qd, lmax_d) = q_method(&obs);
        let aq = quat_to_dcm(qr.q);
        let ad = quat_to_dcm(qd);
        assert!(error_angle(&aq, &ad) < 1e-9, "QUEST and q-method disagree");
        assert!((qr.lambda_max - lmax_d).abs() < 1e-12);
        assert!(qr.newton_iters <= 5);
        // With sum(a_i) = 1, the minimum Wahba loss is 1 - lambda_max.
        assert!((wahba_loss(&obs, &aq) - (1.0 - qr.lambda_max)).abs() < 1e-12);
    }
}

#[test]
fn q_method_solution_is_a_minimum() {
    let mut rng = Rng::new(4);
    let a = quat_to_dcm(random_quat(&mut rng));
    let obs = noisy_set(&mut rng, &a, 4, 0.02);
    let a_opt = quat_to_dcm(q_method(&obs).0);
    let l_opt = wahba_loss(&obs, &a_opt);
    for _ in 0..200 {
        let d = quat_to_dcm(quat_from_axis_angle(rng.unit_vector(), 1e-3));
        assert!(wahba_loss(&obs, &mat_mul(&d, &a_opt)) > l_opt);
    }
}

#[test]
fn quest_survives_180_degree_rotations() {
    let mut rng = Rng::new(5);
    let axes = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 2.0, -1.0], [-3.0, 0.5, 2.0]];
    for &ax in &axes {
        for &ang in &[PI, PI - 1e-6, 3.0] {
            let a = quat_to_dcm(quat_from_axis_angle(ax, ang));
            let obs: Vec<Obs> = (0..3)
                .map(|_| {
                    let r = rng.unit_vector();
                    Obs { b: mat_vec(&a, r), r, sigma: 1e-3 }
                })
                .collect();
            let qr = quest(&obs);
            assert!(error_angle(&quat_to_dcm(qr.q), &a) < 1e-8, "axis {:?} angle {}", ax, ang);
        }
    }
}

#[test]
fn covariance_hand_case() {
    // Two orthogonal unit vectors with equal sigma:
    // F = [(I - x x^T) + (I - y y^T)] / s^2 = diag(1, 1, 2) / s^2
    let s = 0.01;
    let obs = [
        Obs { b: [1.0, 0.0, 0.0], r: [1.0, 0.0, 0.0], sigma: s },
        Obs { b: [0.0, 1.0, 0.0], r: [0.0, 1.0, 0.0], sigma: s },
    ];
    let p = quest_covariance(&obs);
    let want = [s * s, s * s, 0.5 * s * s];
    for i in 0..3 {
        assert!((p[i][i] - want[i]).abs() < 1e-18);
    }
    // TRIAD ignores the second vector's information about the axes the
    // anchor already fixes, so about z it only gets sigma, not sigma/sqrt(2).
    let pt = triad_covariance(obs[0].b, s, obs[1].b, s);
    assert!((pt[2][2] - s * s).abs() < 1e-18);
}

fn monte_carlo_sigma(trials: usize, f: impl Fn(&mut Rng) -> V3) -> V3 {
    let mut rng = Rng::new(6);
    let mut acc = [0.0; 3];
    for _ in 0..trials {
        let e = f(&mut rng);
        for k in 0..3 {
            acc[k] += e[k] * e[k];
        }
    }
    [(acc[0] / trials as f64).sqrt(), (acc[1] / trials as f64).sqrt(), (acc[2] / trials as f64).sqrt()]
}

#[test]
fn monte_carlo_matches_quest_covariance() {
    let sig = [1e-4, 3e-4, 2e-4];
    let body = [normalize([0.1, 0.0, 1.0]), normalize([0.0, 0.2, 1.0]), normalize([-0.15, -0.1, 1.0])];
    let nominal: Vec<Obs> = (0..3).map(|i| Obs { b: body[i], r: body[i], sigma: sig[i] }).collect();
    let p = quest_covariance(&nominal);
    let mc = monte_carlo_sigma(20_000, |rng| {
        let a = quat_to_dcm(random_quat(rng));
        let at = transpose(&a);
        let obs: Vec<Obs> = (0..3)
            .map(|i| Obs { b: measure(rng, body[i], sig[i]), r: mat_vec(&at, body[i]), sigma: sig[i] })
            .collect();
        error_vector(&quat_to_dcm(quest(&obs).q), &a)
    });
    for k in 0..3 {
        let pred = p[k][k].sqrt();
        assert!((mc[k] / pred - 1.0).abs() < 0.03, "axis {}: mc {} pred {}", k, mc[k], pred);
    }
}

#[test]
fn monte_carlo_matches_triad_covariance() {
    let (s1, s2) = (2e-4, 5e-3);
    let b1 = normalize([0.3, -0.2, 1.0]);
    let b2 = normalize([1.0, 0.4, 0.2]);
    let p = triad_covariance(b1, s1, b2, s2);
    let mc = monte_carlo_sigma(20_000, |rng| {
        let a = quat_to_dcm(random_quat(rng));
        let at = transpose(&a);
        let m1 = measure(rng, b1, s1);
        let m2 = measure(rng, b2, s2);
        error_vector(&triad(m1, m2, mat_vec(&at, b1), mat_vec(&at, b2)), &a)
    });
    for k in 0..3 {
        let pred = p[k][k].sqrt();
        assert!((mc[k] / pred - 1.0).abs() < 0.03, "axis {}: mc {} pred {}", k, mc[k], pred);
    }
}
