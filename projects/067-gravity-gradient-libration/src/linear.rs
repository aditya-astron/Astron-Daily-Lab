// Astron Built this AP
// Project 067 | Attitude control | Gravity-gradient stability regions and eccentric-orbit pitch libration
// Author: Aditya (aditya-astron)
// Date: 2026-11-27 IST
//
// Linearised gravity-gradient librations about the LVLH-aligned attitude in a
// circular orbit (time in 1/w0). With the Hughes ratios
//   k1 = (I2 - I3)/I1,  k2 = (I1 - I3)/I2,  k3 = (I2 - I1)/I3
// pitch decouples:
//   theta'' + 3 k2 theta = 0
// and roll/yaw are gyroscopically coupled (signs for this module's frame,
// y = -orbit normal; confirmed against the nonlinear model in the tests):
//   phi'' - (1 - k1) psi' + 4 k1 phi = 0
//   psi'' + (1 - k3) phi' + k3 psi   = 0
// Characteristic equation: s^4 + b s^2 + c = 0 with
//   b = 1 + 3 k1 + k1 k3,  c = 4 k1 k3.
// Both s^2 roots are real and negative (pure oscillation) iff
//   b > 0,  c > 0,  b^2 - 4c > 0.
// Pitch needs k2 > 0, i.e. k1 > k3. The upper-right triangle k1 > k3 > 0 is
// the Lagrange region; the small k1, k3 < 0 lobe is DeBra-Delp, gyroscopic
// stability without a potential minimum.

pub fn k2_from(k1: f64, k3: f64) -> f64 {
    (k1 - k3) / (1.0 - k1 * k3)
}

pub fn pitch_freq(k2: f64) -> Option<f64> {
    if k2 > 0.0 {
        Some((3.0 * k2).sqrt())
    } else {
        None
    }
}

/// Roots of s^4 + b s^2 + c as s^2 values, complex (re, im).
pub fn roll_yaw_s2(k1: f64, k3: f64) -> [(f64, f64); 2] {
    let b = 1.0 + 3.0 * k1 + k1 * k3;
    let c = 4.0 * k1 * k3;
    let d = b * b - 4.0 * c;
    if d >= 0.0 {
        let sq = d.sqrt();
        // stable form: avoid cancellation in the smaller root
        let big = -0.5 * (b + b.signum() * sq);
        let small = if big != 0.0 { c / big } else { 0.0 };
        [(big, 0.0), (small, 0.0)]
    } else {
        let sq = (-d).sqrt();
        [(-0.5 * b, 0.5 * sq), (-0.5 * b, -0.5 * sq)]
    }
}

/// Four characteristic exponents s (re, im) of the roll/yaw system.
pub fn roll_yaw_roots(k1: f64, k3: f64) -> [(f64, f64); 4] {
    let mut out = [(0.0, 0.0); 4];
    for (n, &(re, im)) in roll_yaw_s2(k1, k3).iter().enumerate() {
        // complex sqrt
        let m = (re * re + im * im).sqrt();
        let sr = ((m + re) / 2.0).max(0.0).sqrt();
        let si = ((m - re) / 2.0).max(0.0).sqrt() * if im < 0.0 { -1.0 } else { 1.0 };
        out[2 * n] = (sr, si);
        out[2 * n + 1] = (-sr, -si);
    }
    out
}

/// Roll/yaw natural frequencies (units of w0) if the pair is stable.
pub fn roll_yaw_freqs(k1: f64, k3: f64) -> Option<(f64, f64)> {
    let b = 1.0 + 3.0 * k1 + k1 * k3;
    let c = 4.0 * k1 * k3;
    if b > 0.0 && c > 0.0 && b * b - 4.0 * c > 0.0 {
        let s2 = roll_yaw_s2(k1, k3);
        let (a, b) = ((-s2[0].0).sqrt(), (-s2[1].0).sqrt());
        Some((a.max(b), a.min(b)))
    } else {
        None
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Region {
    Lagrange,
    DeBraDelp,
    Unstable,
}

pub fn classify(k1: f64, k3: f64) -> Region {
    let pitch_ok = k2_from(k1, k3) > 0.0;
    let ry_ok = roll_yaw_freqs(k1, k3).is_some();
    if !(pitch_ok && ry_ok) {
        Region::Unstable
    } else if k1 > 0.0 && k3 > 0.0 {
        Region::Lagrange
    } else {
        Region::DeBraDelp
    }
}

/// State-space roll/yaw: x = (phi, psi, phi', psi').
pub fn roll_yaw_rhs(k1: f64, k3: f64) -> impl Fn(f64, &[f64], &mut [f64]) {
    move |_t, x, dx| {
        dx[0] = x[2];
        dx[1] = x[3];
        dx[2] = (1.0 - k1) * x[3] - 4.0 * k1 * x[0];
        dx[3] = -(1.0 - k3) * x[2] - k3 * x[1];
    }
}
