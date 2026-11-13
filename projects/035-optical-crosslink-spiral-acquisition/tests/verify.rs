// Astron Built this AP
// Project 035 | Laser/optical communication | Beacon spiral-scan acquisition for a LEO-GEO optical relay link
// Author: Aditya (aditya-astron)
// Date: 2026-11-13 IST

use pat_acquisition::reference::*;
use pat_acquisition::*;
use std::f64::consts::{E, PI};

fn rel(a: f64, b: f64) -> f64 { (a - b).abs() / b.abs() }

fn theta_d_ref() -> f64 {
    let l = link();
    l.detection_radius(l.optimal_divergence())
}

#[test]
fn gaussian_tail_reference_values() {
    // tabulated: erfc(1) = 0.157299207, Q^-1(1e-9) = 5.997807, Q^-1(0.01) = 2.326348
    assert!((erfc(1.0) - 0.157_299_207).abs() < 2e-7);
    assert!((erfc(-0.5) - 1.520_499_878).abs() < 2e-7);
    assert!((q_inv(1e-9) - 5.997_807).abs() < 1e-4);
    assert!((q_inv(0.01) - 2.326_348).abs() < 1e-4);
    assert!((q_inv(0.5)).abs() < 1e-6);
}

#[test]
fn threshold_and_received_power_by_hand() {
    let l = link();
    // 8.324 x 4e-14 W/rtHz x sqrt(1 / 20 ms) = 2.354e-12 W
    let snr = 5.997_807 + 2.326_348;
    assert!(rel(l.p_threshold(), snr * 4e-14 * 50f64.sqrt()) < 1e-4);
    // on-axis irradiance 2P/(pi w^2 L^2) times aperture area times optics
    let w = 1e-3;
    let irr = 2.0 * 4.0 / (PI * (w * 45e6) * (w * 45e6));
    let p_rx = irr * PI * 0.135 * 0.135 / 4.0 * 0.45;
    assert!(rel(l.p_rx(w, 0.0), p_rx) < 1e-12);
    // 1/e^2 point of the far field
    assert!(rel(l.p_rx(w, w), p_rx * (-2.0f64).exp()) < 1e-12);
}

#[test]
fn optimal_divergence_puts_e_of_margin_on_axis() {
    let l = link();
    let w_opt = l.optimal_divergence();
    // brute-force scan of theta_d(w)
    let (mut best_w, mut best_th) = (0.0, 0.0);
    for i in 1..200_000 {
        let w = w_opt * 3.0 * i as f64 / 200_000.0;
        let th = l.detection_radius(w);
        if th > best_th { best_th = th; best_w = w; }
    }
    assert!(rel(best_w, w_opt) < 1e-3);
    assert!(rel(l.p_rx(w_opt, 0.0) / l.p_threshold(), E) < 1e-12);
    assert!(rel(l.detection_radius(w_opt), w_opt / 2f64.sqrt()) < 1e-12);
    // the detection edge really is the threshold
    assert!(rel(l.p_rx(w_opt, best_th), l.p_threshold()) < 1e-6);
}

#[test]
fn scan_design_maximises_area_rate() {
    let th = 800e-6;
    let tau = 10e-3;
    for v_max in [1.0, 0.1, 0.05] {
        let d = design_scan(th, tau, v_max);
        // brute force over speed; the widest pitch the dwell allows at each
        let mut best = 0.0f64;
        for i in 1..=100_000 {
            let v = v_max.min(2.0 * th / tau) * i as f64 / 100_000.0;
            let half = 0.5 * v * tau;
            if half >= th { continue; }
            best = best.max(2.0 * (th * th - half * half).sqrt() * v);
        }
        assert!(rel(d.area_rate(), best) < 1e-6, "v_max {v_max}");
        // the worst-case chord equals exactly one dwell
        let chord = 2.0 * (th * th - d.pitch * d.pitch / 4.0).sqrt();
        assert!(rel(chord, d.speed * tau) < 1e-12);
    }
    let free = design_scan(th, tau, 1.0);
    assert!(!free.rate_limited && rel(free.pitch, 2f64.sqrt() * th) < 1e-12);
    assert!(design_scan(th, tau, 0.05).rate_limited);
}

#[test]
fn rayleigh_scan_radius() {
    let r = scan_radius(1.0, 0.999);
    assert!(rel(r, (2.0 * 1000f64.ln()).sqrt()) < 1e-12);
    assert!(rel(1.0 - (-r * r / 2.0).exp(), 0.999) < 1e-12);
}

#[test]
fn spiral_arc_length_and_timing() {
    let sp = Spiral::new(1.2e-3, 0.1, 0.6);
    // trapezoid of |dr/dphi, r| = c sqrt(1 + phi^2)
    let phi_end = 40.0;
    let n = 400_000;
    let mut s = 0.0;
    for i in 0..n {
        let a = phi_end * i as f64 / n as f64;
        let b = phi_end * (i + 1) as f64 / n as f64;
        s += 0.5 * (b - a) * sp.c * ((1.0 + a * a).sqrt() + (1.0 + b * b).sqrt());
    }
    assert!(rel(sp.arc_length(phi_end), s) < 1e-9);
    for phi in [0.01, 1.0, 6.0, 100.0] {
        assert!(rel(sp.phi_at_arc(sp.arc_length(phi)), phi) < 1e-10);
    }
    // beam speed along the track: inner turn, then cruise
    let t1 = sp.arc_length(2.0 * PI) / (0.6 * 0.1);
    for (t, v) in [(0.3 * t1, 0.06), (t1 + 0.5, 0.1), (t1 + 3.0, 0.1)] {
        let h = 1e-6;
        let a = sp.position(t);
        let b = sp.position(t + h);
        assert!(rel((b[0] - a[0]).hypot(b[1] - a[1]) / h, v) < 1e-4, "t {t}");
    }
    assert!(rel(sp.time_to_radius(sp.c * 2.0 * PI), t1) < 1e-12);
}

#[test]
fn windowed_search_matches_brute_force() {
    let th = theta_d_ref();
    let tau = link().tau;
    let d = design_scan(th, tau, V_MAX);
    let track = SpiralTrack::new(&Spiral::new(d.pitch, d.speed, INNER), 6e-3, tau / 50.0);
    let mut rng = Rng::new(7);
    let mut hits = 0;
    for _ in 0..3000 {
        let (gx, gy) = rng.normal2();
        let x = [2.5e-3 * gx, 2.5e-3 * gy];
        let a = track.detect(x, th, tau);
        assert_eq!(a, track.detect_brute(x, th, tau), "target {:?}", x);
        hits += a.is_some() as usize;
    }
    assert!(hits > 2000);
}

#[test]
fn far_field_dwell_tends_to_straight_chord() {
    // A target midway between arms gets its longest dwell from the outer arm,
    // which curves towards it: 2 sqrt(theta_d^2 - d^2) / v stretched by
    // sqrt(r_out / (r_out - d)) for an arc of radius r_out at offset d. With
    // the dwell-limited design the straight chord is exactly tau, so the
    // dwell falls towards tau as 1 + p / (4 r).
    let th = theta_d_ref();
    let tau = link().tau;
    let d = design_scan(th, tau, 10.0);
    let sp = Spiral::new(d.pitch, d.speed, 1.0);
    let track = SpiralTrack::new(&sp, 40e-3, tau / 400.0);
    let mut last = f64::MAX;
    for phi in [40.0, 80.0, 160.0, 220.0] {
        let r = sp.c * phi + 0.5 * d.pitch;
        let r_out = sp.c * (phi + 2.0 * PI);
        let expect = tau * (r_out / (r_out - 0.5 * d.pitch)).sqrt();
        let dwell = track.longest_dwell([r * phi.cos(), r * phi.sin()], th);
        assert!(rel(dwell, expect) < 0.005, "phi {phi}: dwell {dwell} vs {expect}");
        assert!(dwell < last && dwell > 0.99 * tau);
        last = dwell;
    }
    assert!(rel(last, tau) < 0.01);
}

#[test]
fn certified_spiral_covers_the_disk() {
    let th = theta_d_ref();
    let l = link();
    let r_scan = scan_radius(SIGMA, P_UNC);
    let dt = l.tau / 100.0;
    let ideal = design_scan(th, l.tau, V_MAX);
    let cert = certify_pitch(&ideal, INNER, l.tau, r_scan, dt);
    let track = SpiralTrack::new(&Spiral::new(cert.pitch, cert.speed, INNER), r_scan, dt);
    // grids two and four times finer than the one used to certify
    assert!(covers(&track, th, l.tau, r_scan, th / 40.0));
    assert!(covers(&track, th, l.tau, r_scan, th / 80.0));
    let ratio = cert.pitch / ideal.pitch;
    assert!(ratio > 0.95 && ratio < 1.0, "ratio {ratio}");
    // the straight-track pitch flown at constant speed leaves holes
    let naive = SpiralTrack::new(&Spiral::new(ideal.pitch, ideal.speed, 1.0), r_scan, dt);
    assert!(!covers(&naive, th, l.tau, r_scan, th / 20.0));
}

#[test]
fn hex_step_stare_covers_the_disk() {
    let th = theta_d_ref();
    let l = link();
    let r_scan = scan_radius(SIGMA, P_UNC);
    let ss = StepStare::new(th, r_scan, l.tau, T_SETTLE);
    let n = 300;
    for i in 0..=n {
        for j in 0..=n {
            let x = r_scan * (2.0 * i as f64 / n as f64 - 1.0);
            let y = r_scan * (2.0 * j as f64 / n as f64 - 1.0);
            if x.hypot(y) <= r_scan { assert!(ss.detect([x, y]).is_some()); }
        }
    }
    // spot count close to the disk area over the lattice cell area
    let expect = PI * (r_scan + th).powi(2) / ss.cell_area();
    assert!(rel(ss.spots.len() as f64, expect) < 0.06, "{} vs {}", ss.spots.len(), expect);
}

#[test]
fn monte_carlo_matches_models() {
    let l = link();
    let th = theta_d_ref();
    let r_scan = scan_radius(SIGMA, P_UNC);
    let dt = l.tau / 100.0;
    let scan = certify_pitch(&design_scan(th, l.tau, V_MAX), INNER, l.tau, r_scan, dt);
    let spiral = Spiral::new(scan.pitch, scan.speed, INNER);
    let track = SpiralTrack::new(&spiral, r_scan, dt);
    let ss = StepStare::new(th, r_scan, l.tau, T_SETTLE);
    let n = 40_000;

    let mc = monte_carlo(SIGMA, n, 11, |x| track.detect(x, th, l.tau));
    // coverage is at least the designed 0.999 (3 binomial sigmas)
    let sd = (P_UNC * (1.0 - P_UNC) / n as f64).sqrt();
    assert!(mc.p_acq() > P_UNC - 3.0 * sd, "spiral P {}", mc.p_acq());
    let m = spiral_mean_time(&spiral, SIGMA, r_scan + scan.pitch);
    assert!(rel(mc.mean(), m) < 0.04, "spiral mean {} vs {}", mc.mean(), m);

    let mc = monte_carlo(SIGMA, n, 12, |x| ss.detect(x));
    assert!(mc.p_acq() > P_UNC - 3.0 * sd, "stare P {}", mc.p_acq());
    let (_, m) = exponential_model(ss.cell_area() / ss.t_spot, SIGMA, ss.duration());
    assert!(rel(mc.mean(), m) < 0.04, "stare mean {} vs {}", mc.mean(), m);
}

#[test]
fn rng_moments() {
    let mut rng = Rng::new(1);
    let n = 200_000;
    let (mut s1, mut s2) = (0.0, 0.0);
    for _ in 0..n / 2 {
        let (a, b) = rng.normal2();
        s1 += a + b;
        s2 += a * a + b * b;
    }
    assert!((s1 / n as f64).abs() < 0.01);
    assert!((s2 / n as f64 - 1.0).abs() < 0.01);
}
