// Astron Built this AP
// Project 029 | Attitude control | B-dot magnetic detumbling of a 3U CubeSat
// Author: Aditya (aditya-astron)
// Date: 2026-11-11 IST

use bdot_detumble::*;

const D2R: f64 = std::f64::consts::PI / 180.0;

fn cubesat() -> Spacecraft {
    // 4 kg 3U, 10 x 10 x 34 cm, z along the long axis
    Spacecraft { inertia: [0.0405, 0.0405, 0.0067], m_max: 0.2 }
}

fn sso_orbit() -> Orbit {
    Orbit { radius: R_EARTH + 500e3, incl: 97.4 * D2R, raan: 30.0 * D2R, u0: 0.0 }
}

fn reference(k: f64) -> Scenario {
    Scenario {
        sc: cubesat(),
        orbit: sso_orbit(),
        field: Dipole::igrf2020(),
        law: Law::Sampled { k, period: 1.0, quiet: 0.05, noise: 50e-9 },
        w0: [10.0 * D2R, -8.0 * D2R, 6.0 * D2R],
        q0: [0.9, 0.1, -0.3, 0.3],
        dt: 0.05,
        seed: 2029,
    }
}

fn main() {
    let sc = cubesat();
    let orbit = sso_orbit();
    let n = orbit.mean_motion();
    let tp = orbit.period();
    let xi = std::f64::consts::PI - orbit.incl; // fold 97.4 deg to 82.6 deg
    let k_star = avanzini_gain(n, xi, sc.j_min());

    println!("B-dot detumble, 3U CubeSat in a 500 km SSO");
    println!("  inertia {:?} kg m^2, magnetorquers {} A m^2 per axis", sc.inertia, sc.m_max);
    println!("  orbit period {:.1} min, n = {:.4e} rad/s ({:.4} deg/s)", tp / 60.0, n, n / D2R);
    println!("  IGRF-2020 dipole B0 = {:.0} nT, tilt {:.2} deg",
        Dipole::igrf2020().b0 * 1e9,
        (-Dipole::igrf2020().axis_ecef[2]).acos() / D2R);
    println!("  law m = -(k/|B|^2) Bdot, 1 Hz, torquers off 0.1 s around each read, 50 nT noise");
    println!("  Avanzini gain k* = 2n(1+sin xi)Jmin = {:.3e} N m s (xi = {:.1} deg)\n", k_star, xi / D2R);

    let s = reference(k_star);
    let w0 = norm(s.w0);
    let log = simulate(&s, 6.0 * tp, 10.0);
    println!("Reference run, w0 = (10, -8, 6) deg/s, |w0| = {:.2} deg/s", w0 / D2R);
    println!("  {:>7} {:>10} {:>10} {:>10} {:>10} {:>11}", "orbits", "wx deg/s", "wy deg/s", "wz deg/s", "|w| deg/s", "E (uJ)");
    let marks = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0];
    for m in marks {
        let idx = ((m * tp / 10.0).round() as usize).min(log.len() - 1);
        let p = &log[idx];
        println!("  {:>7.2} {:>10.3} {:>10.3} {:>10.3} {:>10.3} {:>11.3}",
            p.t / tp, p.w[0] / D2R, p.w[1] / D2R, p.w[2] / D2R, norm(p.w) / D2R, p.energy * 1e6);
    }
    let sat_frac = log.iter().filter(|p| p.dipole > 0.999 * sc.m_max).count() as f64 / log.len() as f64;
    println!("  samples at saturation: {:.1} %", 100.0 * sat_frac);
    for th in [1.0, 0.5] {
        match settle_time(&log, th * D2R) {
            Some(t) => println!("  |w| < {:.1} deg/s for good after {:.2} orbits ({:.0} min)", th, t / tp, t / 60.0),
            None => println!("  |w| < {:.1} deg/s not reached", th),
        }
    }
    let tail: Vec<f64> = log.iter().filter(|p| p.t > 5.0 * tp).map(|p| norm(p.w) / D2R).collect();
    let tail_mean = tail.iter().sum::<f64>() / tail.len() as f64;
    println!("  residual |w| over orbit 6: mean {:.3} deg/s (2n = {:.3} deg/s)\n", tail_mean, 2.0 * n / D2R);

    // k*_max: the same formula evaluated with the largest inertia.
    let jr = sc.j_min() / sc.inertia.iter().cloned().fold(0.0, f64::max);
    println!("Gain sweep (same initial state, time to |w| < 0.5 deg/s for good)");
    println!("  {:>8} {:>9} {:>12} {:>10}", "k/k*", "k/k*_max", "k (N m s)", "orbits");
    for f in [0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0] {
        let log = simulate(&reference(f * k_star), 8.0 * tp, 10.0);
        let t = settle_time(&log, 0.5 * D2R).map(|t| format!("{:.2}", t / tp)).unwrap_or("> 8".into());
        println!("  {:>8.3} {:>9.3} {:>12.3e} {:>10}", f, f * jr, f * k_star, t);
    }

    println!("\nSame sweep with the ideal law (no saturation, exact rate):");
    println!("  {:>8} {:>10} {:>14}", "k/k*", "orbits", "peak m (A m^2)");
    for f in [0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0] {
        let mut s = reference(f * k_star);
        s.law = Law::Ideal { k: f * k_star };
        let log = simulate(&s, 8.0 * tp, 10.0);
        let t = settle_time(&log, 0.5 * D2R).map(|t| format!("{:.2}", t / tp)).unwrap_or("> 8".into());
        let peak = log.iter().fold(0.0f64, |a, p| a.max(p.dipole));
        println!("  {:>8.3} {:>10} {:>14.2}", f, t, peak);
    }

    println!("\nAveraged theory, aligned dipole, isotropic J: slowest rate (k/J)(1 - lambda_max)");
    println!("  {:>8} {:>11} {:>15} {:>20}", "xi deg", "lambda_max", "rate/(k/J)", "rate at k* / n");
    for xi_deg in [0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0] {
        let ev = sym_eigenvalues(averaged_bbt(xi_deg * D2R, 720));
        let r = 1.0 - ev[2];
        println!("  {:>8.0} {:>11.4} {:>15.4} {:>20.4}", xi_deg, ev[2], r,
            r * 2.0 * (1.0 + (xi_deg * D2R).sin()));
    }
}
