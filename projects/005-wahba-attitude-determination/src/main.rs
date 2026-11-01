// Astron Built this AP
// Project 005 | Attitude determination | Wahba's problem: TRIAD, Davenport q-method and QUEST
// Author: Aditya (aditya-astron)
// Date: 2026-11-01 IST
//
// Monte Carlo comparison of three static attitude solvers on two sensor
// suites, checked against the analytic covariances.

use std::time::Instant;
use wahba::*;

const ARCSEC: f64 = std::f64::consts::PI / (180.0 * 3600.0);
const DEG: f64 = std::f64::consts::PI / 180.0;
const TRIALS: usize = 20_000;

/// Accumulates the sample covariance of body-frame error vectors.
struct Stats {
    sum: M3,
    n: usize,
}

impl Stats {
    fn new() -> Self {
        Stats { sum: [[0.0; 3]; 3], n: 0 }
    }
    fn push(&mut self, e: V3) {
        self.sum = mat_add(&self.sum, &outer(e, e));
        self.n += 1;
    }
    fn sigma(&self) -> V3 {
        let c = mat_scale(&self.sum, 1.0 / self.n as f64);
        [c[0][0].sqrt(), c[1][1].sqrt(), c[2][2].sqrt()]
    }
}

fn diag_sqrt(p: &M3) -> V3 {
    [p[0][0].sqrt(), p[1][1].sqrt(), p[2][2].sqrt()]
}

fn row(label: &str, s: V3, unit: f64) {
    println!(
        "  {:<28} {:>9.3} {:>9.3} {:>9.3}",
        label,
        s[0] / unit,
        s[1] / unit,
        s[2] / unit
    );
}

fn random_attitude(rng: &mut Rng) -> M3 {
    quat_to_dcm(quat_normalize([rng.gauss(), rng.gauss(), rng.gauss(), rng.gauss()]))
}

fn star_tracker(rng: &mut Rng) {
    // Six stars inside an 8 deg half-angle field of view around body +z.
    // Cross-boresight accuracy of 3 arcsec per star is typical for a
    // modern star tracker; roll about the boresight is only observable
    // through the small lever arm of the star spread, so it comes out much
    // worse.
    let sigma = 3.0 * ARCSEC;
    let half_fov = 8.0 * DEG;
    let stars: Vec<V3> = (0..6)
        .map(|_| {
            let rho = half_fov * rng.uniform().sqrt();
            let phi = 2.0 * std::f64::consts::PI * rng.uniform();
            [rho.sin() * phi.cos(), rho.sin() * phi.sin(), rho.cos()]
        })
        .collect();
    let nominal: Vec<Obs> = stars.iter().map(|&b| Obs { b, r: b, sigma }).collect();
    let p = quest_covariance(&nominal);

    let mut st_q = Stats::new();
    let mut st_d = Stats::new();
    let mut max_disagree: f64 = 0.0;
    let mut iters_max = 0;
    for _ in 0..TRIALS {
        let a_true = random_attitude(rng);
        let at = transpose(&a_true);
        let obs: Vec<Obs> = stars
            .iter()
            .map(|&b| Obs { b: measure(rng, b, sigma), r: mat_vec(&at, b), sigma })
            .collect();
        let qr = quest(&obs);
        let (qd, _) = q_method(&obs);
        let aq = quat_to_dcm(qr.q);
        let ad = quat_to_dcm(qd);
        st_q.push(error_vector(&aq, &a_true));
        st_d.push(error_vector(&ad, &a_true));
        max_disagree = max_disagree.max(error_angle(&aq, &ad));
        iters_max = iters_max.max(qr.newton_iters);
    }

    println!("Case 1: star tracker, 6 stars, 8 deg half-FOV, 3 arcsec/star, {} trials", TRIALS);
    println!("  1-sigma error, arcsec         {:>9} {:>9} {:>9}", "x", "y", "z(roll)");
    row("predicted (QUEST covariance)", diag_sqrt(&p), ARCSEC);
    row("QUEST", st_q.sigma(), ARCSEC);
    row("q-method", st_d.sigma(), ARCSEC);
    println!(
        "  QUEST vs q-method max difference: {:.2e} arcsec, Newton iterations <= {}",
        max_disagree / ARCSEC,
        iters_max
    );
    println!();
}

fn sun_mag(rng: &mut Rng) {
    // Fine sun sensor (0.05 deg) and magnetometer (1 deg after field-model
    // error), 60 deg apart in the body frame.
    let s_sun = 0.05 * DEG;
    let s_mag = 1.0 * DEG;
    let b_sun = [0.0, 0.0, 1.0];
    let b_mag = [(60.0 * DEG).sin(), 0.0, (60.0 * DEG).cos()];

    let nominal = [
        Obs { b: b_sun, r: b_sun, sigma: s_sun },
        Obs { b: b_mag, r: b_mag, sigma: s_mag },
    ];
    let p_quest = quest_covariance(&nominal);
    let p_t_sun = triad_covariance(b_sun, s_sun, b_mag, s_mag);
    let p_t_mag = triad_covariance(b_mag, s_mag, b_sun, s_sun);

    let mut st_ts = Stats::new();
    let mut st_tm = Stats::new();
    let mut st_q = Stats::new();
    for _ in 0..TRIALS {
        let a_true = random_attitude(rng);
        let at = transpose(&a_true);
        let r_sun = mat_vec(&at, b_sun);
        let r_mag = mat_vec(&at, b_mag);
        let m_sun = measure(rng, b_sun, s_sun);
        let m_mag = measure(rng, b_mag, s_mag);
        st_ts.push(error_vector(&triad(m_sun, m_mag, r_sun, r_mag), &a_true));
        st_tm.push(error_vector(&triad(m_mag, m_sun, r_mag, r_sun), &a_true));
        let obs = [
            Obs { b: m_sun, r: r_sun, sigma: s_sun },
            Obs { b: m_mag, r: r_mag, sigma: s_mag },
        ];
        st_q.push(error_vector(&quat_to_dcm(quest(&obs).q), &a_true));
    }

    println!("Case 2: sun sensor 0.05 deg + magnetometer 1 deg, 60 deg apart, {} trials", TRIALS);
    println!("  1-sigma error, deg            {:>9} {:>9} {:>9}", "x", "y", "z(sun)");
    row("TRIAD sun-first  predicted", diag_sqrt(&p_t_sun), DEG);
    row("TRIAD sun-first  simulated", st_ts.sigma(), DEG);
    row("TRIAD mag-first  predicted", diag_sqrt(&p_t_mag), DEG);
    row("TRIAD mag-first  simulated", st_tm.sigma(), DEG);
    row("QUEST            predicted", diag_sqrt(&p_quest), DEG);
    row("QUEST            simulated", st_q.sigma(), DEG);
    println!();
}

fn singular_case() {
    // 180 deg rotation: plain QUEST's gamma goes to zero here.
    let q_true = quat_from_axis_angle([1.0, 2.0, -1.0], std::f64::consts::PI);
    let a_true = quat_to_dcm(q_true);
    let refs = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], normalize([1.0, 1.0, 1.0])];
    let obs: Vec<Obs> = refs
        .iter()
        .map(|&r| Obs { b: mat_vec(&a_true, r), r, sigma: 1e-4 })
        .collect();
    let qr = quest(&obs);
    println!("Case 3: exact 180 deg rotation (QUEST singularity)");
    println!(
        "  sequential rotation about axis {:?}, attitude error {:.2e} rad",
        qr.rotated_axis.map(|a| ["x", "y", "z"][a]),
        error_angle(&quat_to_dcm(qr.q), &a_true)
    );
    println!();
}

fn timing(rng: &mut Rng) {
    let n = 200_000;
    let sets: Vec<Vec<Obs>> = (0..256)
        .map(|_| {
            let a = random_attitude(rng);
            (0..6)
                .map(|_| {
                    let r = rng.unit_vector();
                    Obs { b: measure(rng, mat_vec(&a, r), 1e-4), r, sigma: 1e-4 }
                })
                .collect()
        })
        .collect();
    let mut sink = 0.0;
    let t0 = Instant::now();
    for i in 0..n {
        sink += quest(&sets[i & 255]).q[3];
    }
    let tq = t0.elapsed().as_secs_f64() / n as f64;
    let t0 = Instant::now();
    for i in 0..n {
        sink += q_method(&sets[i & 255]).0[3];
    }
    let td = t0.elapsed().as_secs_f64() / n as f64;
    println!("Case 4: cost per solve, 6 observations ({} solves, checksum {:.3})", n, sink);
    println!("  QUEST    {:.2} us", tq * 1e6);
    println!("  q-method {:.2} us  (Jacobi eigen-solve)", td * 1e6);
}

fn main() {
    let mut rng = Rng::new(0x5EED_A57A);
    star_tracker(&mut rng);
    sun_mag(&mut rng);
    singular_case();
    timing(&mut rng);
}
