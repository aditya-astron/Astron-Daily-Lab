// Astron Built this AP
// Project 035 | Laser/optical communication | Beacon spiral-scan acquisition for a LEO-GEO optical relay link
// Author: Aditya (aditya-astron)
// Date: 2026-11-13 IST

//! Open-loop acquisition for an optical inter-satellite link.
//!
//! Terminal A (GEO relay) sweeps a Gaussian beacon over its pointing
//! uncertainty cone. Terminal B (LEO user) stares with a wide-FOV acquisition
//! sensor and declares detection once the received beacon power has stayed
//! above threshold for one integration time `tau`. Everything is in the
//! small-angle plane around A's best-estimate line of sight (radians).
//!
//! Design chain:
//!   sensor NEP, tau, Pfa, Pd      -> threshold power P_th
//!   beacon power, range, aperture -> detection radius theta_d(w)
//!   theta_d, tau, slew limit      -> spiral pitch p and scan speed v
//!   pointing sigma, P_unc         -> scan radius R, scan time, E[T_acq]

use std::f64::consts::{E, PI};

// ---------------------------------------------------------------- statistics

/// Complementary error function, Chebyshev fit from Numerical Recipes
/// (fractional error < 1.2e-7 everywhere), good enough for threshold design.
pub fn erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.5 * z);
    let ans = t
        * (-z * z - 1.26551223
            + t * (1.00002368
                + t * (0.37409196
                    + t * (0.09678418
                        + t * (-0.18628806
                            + t * (0.27886807
                                + t * (-1.13520398
                                    + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))))
            .exp();
    if x >= 0.0 { ans } else { 2.0 - ans }
}

/// Gaussian tail Q(x) = P(N(0,1) > x).
pub fn q_func(x: f64) -> f64 {
    0.5 * erfc(x / 2f64.sqrt())
}

/// Inverse of Q by bisection; Q is monotone so this is robust.
pub fn q_inv(p: f64) -> f64 {
    let (mut lo, mut hi) = (-40.0, 40.0);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if q_func(mid) > p { lo = mid } else { hi = mid }
    }
    0.5 * (lo + hi)
}

/// SplitMix64: tiny, well-mixed, reproducible. No external crates needed.
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self { Rng(seed) }
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    /// Uniform in (0, 1).
    pub fn uniform(&mut self) -> f64 {
        ((self.next_u64() >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    }
    /// Pair of independent N(0,1) samples (Box-Muller).
    pub fn normal2(&mut self) -> (f64, f64) {
        let r = (-2.0 * self.uniform().ln()).sqrt();
        let a = 2.0 * PI * self.uniform();
        (r * a.cos(), r * a.sin())
    }
}

// ---------------------------------------------------------------- link model

#[derive(Clone, Copy, Debug)]
pub struct Link {
    pub p_beacon: f64,  // W, beacon optical power out of the telescope
    pub range: f64,     // m
    pub d_rx: f64,      // m, acquisition sensor aperture diameter
    pub eta_rx: f64,    // receive optics x filter x fill factor
    pub nep: f64,       // W/sqrt(Hz), detector noise-equivalent power
    pub tau: f64,       // s, integration (dwell) time
    pub pfa: f64,       // false-alarm probability per decision
    pub pd: f64,        // detection probability per decision
}

impl Link {
    /// Required SNR for a threshold test on Gaussian noise:
    /// threshold sits Q^-1(Pfa) sigmas above zero and the signal must sit
    /// Q^-1(1-Pd) sigmas above the threshold.
    pub fn required_snr(&self) -> f64 {
        q_inv(self.pfa) + q_inv(1.0 - self.pd)
    }

    /// Threshold power: SNR x NEP x sqrt(B), B = 1/(2 tau) for an
    /// integrate-and-dump of length tau.
    pub fn p_threshold(&self) -> f64 {
        self.required_snr() * self.nep * (1.0 / (2.0 * self.tau)).sqrt()
    }

    /// Collected power constant: P_rx(theta) = (K / w^2) exp(-2 theta^2 / w^2)
    /// for a Gaussian far field with 1/e^2 half-angle divergence w, since the
    /// on-axis irradiance is 2 P / (pi w^2 L^2).
    pub fn k_rx(&self) -> f64 {
        let area = PI * self.d_rx * self.d_rx / 4.0;
        2.0 * self.p_beacon * area * self.eta_rx / (PI * self.range * self.range)
    }

    pub fn p_rx(&self, w: f64, theta: f64) -> f64 {
        self.k_rx() / (w * w) * (-2.0 * theta * theta / (w * w)).exp()
    }

    /// Angle out to which P_rx >= P_th (0 if the link closes nowhere).
    pub fn detection_radius(&self, w: f64) -> f64 {
        let margin = self.k_rx() / (w * w * self.p_threshold());
        if margin <= 1.0 { 0.0 } else { w * (0.5 * margin.ln()).sqrt() }
    }

    /// theta_d^2 = (w^2/2) ln(C/w^2), C = K/P_th. d/dw = 0 gives
    /// ln(C/w^2) = 1: the best beacon puts exactly e (4.34 dB) of margin on
    /// axis, and then theta_d = w / sqrt(2).
    pub fn optimal_divergence(&self) -> f64 {
        (self.k_rx() / (E * self.p_threshold())).sqrt()
    }
}

// ---------------------------------------------------------------- scan design

#[derive(Clone, Copy, Debug)]
pub struct ScanDesign {
    pub theta_d: f64, // rad, detection radius
    pub pitch: f64,   // rad, spiral arm spacing
    pub speed: f64,   // rad/s, beam speed along the spiral
    pub rate_limited: bool,
}

impl ScanDesign {
    /// Swept-area rate, rad^2/s.
    pub fn area_rate(&self) -> f64 { self.pitch * self.speed }
}

/// Pitch and speed that maximise swept area per second.
///
/// A target midway between arms sits p/2 off the beam track and is lit for a
/// chord 2 sqrt(theta_d^2 - p^2/4) of track. Dwell needs chord >= v tau, so the
/// area rate is p v <= (2/tau) p sqrt(theta_d^2 - p^2/4), maximised at
/// p = sqrt(2) theta_d, v = sqrt(2) theta_d / tau. If the scanner cannot reach
/// that speed, run at v_max and open the pitch to the widest the dwell allows.
pub fn design_scan(theta_d: f64, tau: f64, v_max: f64) -> ScanDesign {
    let v_free = 2f64.sqrt() * theta_d / tau;
    if v_free <= v_max {
        ScanDesign { theta_d, pitch: 2f64.sqrt() * theta_d, speed: v_free, rate_limited: false }
    } else {
        let half_chord = 0.5 * v_max * tau;
        let pitch = 2.0 * (theta_d * theta_d - half_chord * half_chord).sqrt();
        ScanDesign { theta_d, pitch, speed: v_max, rate_limited: true }
    }
}

/// Radius holding the target with probability p_unc when its offset is a
/// circular Gaussian with per-axis sigma (Rayleigh radius).
pub fn scan_radius(sigma: f64, p_unc: f64) -> f64 {
    sigma * (-2.0 * (1.0 - p_unc).ln()).sqrt()
}

// ---------------------------------------------------------------- spiral

/// Archimedean spiral r = c phi, c = p / (2 pi). The beam runs at `speed`
/// along the track except on the first turn, which is flown at
/// `inner * speed`: there the arm radius is comparable to theta_d, the track
/// curls inside the detection circle and the straight-track dwell argument
/// no longer holds.
pub struct Spiral {
    pub c: f64,
    pub speed: f64,
    pub inner: f64,
}

impl Spiral {
    pub fn new(pitch: f64, speed: f64, inner: f64) -> Self {
        Spiral { c: pitch / (2.0 * PI), speed, inner }
    }

    /// Exact arc length from the centre to angle phi.
    pub fn arc_length(&self, phi: f64) -> f64 {
        0.5 * self.c * (phi * (1.0 + phi * phi).sqrt() + phi.asinh())
    }

    /// Spiral angle reached after travelling s along the track (Newton on the
    /// exact arc length, started from the large-phi asymptote s ~ c phi^2 / 2).
    pub fn phi_at_arc(&self, s: f64) -> f64 {
        if s <= 0.0 { return 0.0; }
        let mut phi = (2.0 * s / self.c).sqrt();
        for _ in 0..50 {
            let f = self.arc_length(phi) - s;
            let step = f / (self.c * (1.0 + phi * phi).sqrt());
            phi -= step;
            if step.abs() < 1e-14 * (1.0 + phi) { break; }
        }
        phi
    }

    fn first_turn(&self) -> (f64, f64) {
        let s1 = self.arc_length(2.0 * PI);
        (s1, s1 / (self.inner * self.speed))
    }

    /// Arc length flown by time t.
    pub fn arc_at_time(&self, t: f64) -> f64 {
        let (s1, t1) = self.first_turn();
        if t <= t1 { self.inner * self.speed * t } else { s1 + self.speed * (t - t1) }
    }

    /// Time for the beam to reach radius r.
    pub fn time_to_radius(&self, r: f64) -> f64 {
        let s = self.arc_length(r / self.c);
        let (s1, t1) = self.first_turn();
        if s <= s1 { s / (self.inner * self.speed) } else { t1 + (s - s1) / self.speed }
    }

    pub fn position(&self, t: f64) -> [f64; 2] {
        let phi = self.phi_at_arc(self.arc_at_time(t));
        let r = self.c * phi;
        [r * phi.cos(), r * phi.sin()]
    }
}

/// Beam track sampled at a fixed step, kept once and reused for every target.
pub struct SpiralTrack {
    pub dt: f64,
    pub pos: Vec<[f64; 2]>,
    pub rad: Vec<f64>, // non-decreasing
    pub phi: Vec<f64>, // unwrapped spiral angle, increasing
}

impl SpiralTrack {
    /// Fly out until the beam centre is a full pitch beyond R. Points just
    /// inside R then have an arm curving towards them on the outside as well;
    /// stopping at R + p/2 leaves them on the convex side of the last arm,
    /// where the dwell is shortest.
    pub fn new(spiral: &Spiral, r_scan: f64, dt: f64) -> Self {
        let pitch = 2.0 * PI * spiral.c;
        let t_end = spiral.time_to_radius(r_scan + pitch);
        let n = (t_end / dt).ceil() as usize + 1;
        let mut pos = Vec::with_capacity(n);
        let mut rad = Vec::with_capacity(n);
        let mut phi = Vec::with_capacity(n);
        for k in 0..n {
            let f = spiral.phi_at_arc(spiral.arc_at_time(k as f64 * dt));
            let r = spiral.c * f;
            pos.push([r * f.cos(), r * f.sin()]);
            rad.push(r);
            phi.push(f);
        }
        SpiralTrack { dt, pos, rad, phi }
    }

    pub fn duration(&self) -> f64 { (self.pos.len() - 1) as f64 * self.dt }

    /// Sample index ranges, in time order, that can possibly light the
    /// target: |r - rho| <= theta_d, and for rho > theta_d also an angular
    /// window of asin(theta_d / rho) around the target's bearing on each arm.
    /// A beam outside these windows is outside the detection circle, so dwell
    /// runs never span two windows.
    pub fn windows(&self, target: [f64; 2], theta_d: f64) -> Vec<(usize, usize)> {
        let rho = target[0].hypot(target[1]);
        let lo = self.rad.partition_point(|&r| r < rho - theta_d);
        let hi = self.rad.partition_point(|&r| r <= rho + theta_d);
        if lo >= hi { return vec![]; }
        if rho <= 2.0 * theta_d { return vec![(lo, hi)]; }
        let half = (theta_d / rho).asin() * (1.0 + 1e-9);
        let bearing = target[1].atan2(target[0]).rem_euclid(2.0 * PI);
        let (f_lo, f_hi) = (self.phi[lo], self.phi[hi - 1]);
        let m0 = ((f_lo - half - bearing) / (2.0 * PI)).floor() as i64;
        let m1 = ((f_hi + half - bearing) / (2.0 * PI)).ceil() as i64;
        let mut out = Vec::new();
        for m in m0.max(0)..=m1 {
            let fc = bearing + 2.0 * PI * m as f64;
            let a = self.phi.partition_point(|&f| f < fc - half).max(lo);
            let b = self.phi.partition_point(|&f| f <= fc + half).min(hi);
            if a < b { out.push((a, b)); }
        }
        out
    }

    fn inside(&self, k: usize, target: [f64; 2], d2: f64) -> bool {
        let dx = self.pos[k][0] - target[0];
        let dy = self.pos[k][1] - target[1];
        dx * dx + dy * dy <= d2
    }

    /// Longest continuous time the target spends inside the beam's
    /// detection radius during the whole scan.
    pub fn longest_dwell(&self, target: [f64; 2], theta_d: f64) -> f64 {
        let d2 = theta_d * theta_d;
        let mut best = 0usize;
        for (a, b) in self.windows(target, theta_d) {
            let mut run = 0usize;
            for k in a..b {
                if self.inside(k, target, d2) { run += 1; best = best.max(run); } else { run = 0; }
            }
        }
        best.saturating_sub(1) as f64 * self.dt
    }

    /// Time at which B declares detection, or None if the scan misses it.
    /// Detection = beam within theta_d of the target continuously for tau.
    pub fn detect(&self, target: [f64; 2], theta_d: f64, tau: f64) -> Option<f64> {
        self.detect_in(&self.windows(target, theta_d), target, theta_d, tau)
    }

    /// Same decision scanning every sample of the radial band; reference
    /// for the windowed search.
    pub fn detect_brute(&self, target: [f64; 2], theta_d: f64, tau: f64) -> Option<f64> {
        self.detect_in(&[(0, self.pos.len())], target, theta_d, tau)
    }

    fn detect_in(&self, ranges: &[(usize, usize)], target: [f64; 2], theta_d: f64, tau: f64) -> Option<f64> {
        let need = (tau / self.dt).ceil() as usize;
        let d2 = theta_d * theta_d;
        for &(a, b) in ranges {
            let mut run = 0usize;
            for k in a..b {
                if self.inside(k, target, d2) {
                    run += 1;
                    // run samples span (run - 1) dt of continuous illumination
                    if run > need { return Some(k as f64 * self.dt); }
                } else {
                    run = 0;
                }
            }
        }
        None
    }
}

/// Every point of a dense polar grid inside r_scan is detected.
pub fn covers(track: &SpiralTrack, theta_d: f64, tau: f64, r_scan: f64, step: f64) -> bool {
    let nr = (r_scan / step).ceil() as usize;
    for i in 0..=nr {
        let r = r_scan * i as f64 / nr as f64;
        let na = ((2.0 * PI * r / step).ceil() as usize).max(1);
        for j in 0..na {
            let a = 2.0 * PI * j as f64 / na as f64;
            if track.detect([r * a.cos(), r * a.sin()], theta_d, tau).is_none() {
                return false;
            }
        }
    }
    true
}

/// The straight-track pitch is exact only where the arms are nearly straight
/// and the target midway between two of them sees the same chord from both.
/// Near the centre and on the convex side of the last arm the lit arc is
/// shorter, and a sampled dwell adds a quantisation loss. Keep the speed and
/// bisect the largest pitch that still covers a dense grid of the whole disk.
/// That pitch sits on the edge by construction, so back off 1 % to cover the
/// points between grid nodes.
pub fn certify_pitch(design: &ScanDesign, inner: f64, tau: f64, r_scan: f64, dt: f64) -> ScanDesign {
    let step = design.theta_d / 20.0;
    let ok = |p: f64| {
        let track = SpiralTrack::new(&Spiral::new(p, design.speed, inner), r_scan, dt);
        covers(&track, design.theta_d, tau, r_scan, step)
    };
    const BACKOFF: f64 = 0.99;
    let (mut lo, mut hi) = (0.3 * design.pitch, design.pitch);
    if ok(hi) { return ScanDesign { pitch: BACKOFF * hi, ..*design }; }
    assert!(ok(lo), "no covering pitch found");
    while hi - lo > 1e-3 * design.pitch {
        let mid = 0.5 * (lo + hi);
        if ok(mid) { lo = mid } else { hi = mid }
    }
    ScanDesign { pitch: BACKOFF * lo, ..*design }
}

// ---------------------------------------------------------------- step-stare

/// Hexagonal step-and-stare: spots of radius theta_d on a triangular lattice
/// with spacing sqrt(3) theta_d tile the plane with no gaps. Each spot costs
/// one settle time plus one dwell; spots are visited ring by ring outward.
pub struct StepStare {
    pub spots: Vec<[f64; 2]>,
    pub t_spot: f64,
    pub theta_d: f64,
}

impl StepStare {
    pub fn new(theta_d: f64, r_scan: f64, tau: f64, t_settle: f64) -> Self {
        let a = 3f64.sqrt() * theta_d;
        let n = (r_scan / a).ceil() as i64 + 2;
        let mut spots = Vec::new();
        for j in -n..=n {
            for i in -2 * n..=2 * n {
                let x = a * (i as f64 + 0.5 * j as f64);
                let y = a * 0.5 * 3f64.sqrt() * j as f64;
                // centres beyond R + theta_d cannot light anything inside R
                if x.hypot(y) <= r_scan + theta_d { spots.push([x, y]); }
            }
        }
        // lattice rings first, then angle: keeps the pattern growing outward
        spots.sort_by(|p, q| {
            let rp = (p[0].hypot(p[1]) / a * 1e6).round();
            let rq = (q[0].hypot(q[1]) / a * 1e6).round();
            rp.partial_cmp(&rq).unwrap().then(p[1].atan2(p[0]).partial_cmp(&q[1].atan2(q[0])).unwrap())
        });
        StepStare { spots, t_spot: t_settle + tau, theta_d }
    }

    /// Lattice cell area per spot.
    pub fn cell_area(&self) -> f64 { 1.5 * 3f64.sqrt() * self.theta_d * self.theta_d }

    pub fn duration(&self) -> f64 { self.spots.len() as f64 * self.t_spot }

    pub fn detect(&self, target: [f64; 2]) -> Option<f64> {
        let d2 = self.theta_d * self.theta_d;
        self.spots.iter().position(|s| {
            let dx = s[0] - target[0];
            let dy = s[1] - target[1];
            dx * dx + dy * dy <= d2
        })
        .map(|k| (k + 1) as f64 * self.t_spot)
    }
}

// ---------------------------------------------------------------- statistics of T_acq

/// Covering a disk at swept-area rate A_dot reaches radius rho after
/// pi rho^2 / A_dot. With rho Rayleigh(sigma), rho^2 is exponential, so T_acq
/// is exponential with rate lambda = A_dot / (2 pi sigma^2), truncated at the
/// end of the scan. Returns (P(acquired by t_end), E[T | acquired]).
pub fn exponential_model(area_rate: f64, sigma: f64, t_end: f64) -> (f64, f64) {
    let lam = area_rate / (2.0 * PI * sigma * sigma);
    let p = 1.0 - (-lam * t_end).exp();
    let mean = 1.0 / lam - t_end * (-lam * t_end).exp() / p;
    (p, mean)
}

#[derive(Debug, Default, Clone)]
pub struct McResult {
    pub n: usize,
    pub hits: usize,
    pub times: Vec<f64>, // acquisition times of the hits, sorted
}

impl McResult {
    pub fn p_acq(&self) -> f64 { self.hits as f64 / self.n as f64 }
    pub fn mean(&self) -> f64 { self.times.iter().sum::<f64>() / self.times.len() as f64 }
    pub fn quantile(&self, q: f64) -> f64 {
        let k = ((self.times.len() as f64 - 1.0) * q).round() as usize;
        self.times[k]
    }
}

/// Draw true target offsets from the pointing error distribution and run a
/// detector over each.
pub fn monte_carlo<F: Fn([f64; 2]) -> Option<f64>>(sigma: f64, n: usize, seed: u64, detect: F) -> McResult {
    let mut rng = Rng::new(seed);
    let mut times = Vec::with_capacity(n);
    for _ in 0..n {
        let (gx, gy) = rng.normal2();
        if let Some(t) = detect([sigma * gx, sigma * gy]) { times.push(t); }
    }
    times.sort_by(|a, b| a.partial_cmp(b).unwrap());
    McResult { n, hits: times.len(), times }
}

/// Semi-analytic spiral acquisition time: the beam reaches radius rho at
/// T(rho) (exact arc length, slow first turn included); average over the
/// Rayleigh radius distribution truncated at r_max (composite Simpson).
pub fn spiral_mean_time(spiral: &Spiral, sigma: f64, r_max: f64) -> f64 {
    let n = 4000;
    let h = r_max / n as f64;
    let pdf = |r: f64| r / (sigma * sigma) * (-r * r / (2.0 * sigma * sigma)).exp();
    let (mut num, mut den) = (0.0, 0.0);
    for i in 0..=n {
        let r = i as f64 * h;
        let w = if i == 0 || i == n { 1.0 } else if i % 2 == 1 { 4.0 } else { 2.0 };
        num += w * pdf(r) * spiral.time_to_radius(r);
        den += w * pdf(r);
    }
    num / den
}

// ---------------------------------------------------------------- reference case

/// GEO data-relay terminal acquiring a LEO user near the edge of its coverage.
pub mod reference {
    use super::Link;

    pub const SIGMA: f64 = 2.5e-3; // rad per axis, GEO open-loop pointing error
    pub const P_UNC: f64 = 0.999; // share of the error distribution to cover
    pub const V_MAX: f64 = 0.10; // rad/s, coarse pointing slew limit
    pub const T_SETTLE: f64 = 4e-3; // s, step-and-stare FSM step and settle
    pub const INNER: f64 = 0.6; // first spiral turn flown at 60 % speed

    pub fn link() -> Link {
        Link {
            p_beacon: 4.0,
            range: 45_000e3,
            d_rx: 0.135,
            eta_rx: 0.45,
            nep: 4e-14,
            tau: 10e-3,
            pfa: 1e-9,
            pd: 0.99,
        }
    }
}
