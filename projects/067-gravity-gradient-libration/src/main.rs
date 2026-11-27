// Astron Built this AP
// Project 067 | Attitude control | Gravity-gradient stability regions and eccentric-orbit pitch libration
// Author: Aditya (aditya-astron)
// Date: 2026-11-27 IST
//
// Design study: a 50 kg microsat stabilised by a nadir boom in a 600 x 880 km
// orbit (e = 0.0197). How long must the boom be to keep the eccentricity-
// forced pitch libration under 2 deg, and does the chosen length sit clear
// of the parametric tongues and stay stable in all three axes?

use gg_libration::linear::{self, Region};
use gg_libration::ode;
use gg_libration::pitch;
use gg_libration::rigid::Body;
use std::f64::consts::PI;

const RE: f64 = 6378.137;

struct Craft {
    bus_mass: f64,
    bus_i: [f64; 3],  // roll, pitch, yaw about bus CoM, kg m^2
    attach: f64,      // boom root offset from bus CoM along yaw axis, m
    rod_density: f64, // kg/m
    tip_mass: f64,
}

impl Craft {
    /// Composite principal inertias for a boom of length L on the yaw axis.
    /// The boom is thin, so it adds only to roll and pitch.
    fn inertia(&self, l: f64) -> [f64; 3] {
        let rod_m = self.rod_density * l;
        let z_rod = self.attach + 0.5 * l;
        let z_tip = self.attach + l;
        let m = self.bus_mass + rod_m + self.tip_mass;
        let zc = (rod_m * z_rod + self.tip_mass * z_tip) / m;
        let add = self.bus_mass * zc * zc
            + rod_m * ((z_rod - zc).powi(2) + l * l / 12.0)
            + self.tip_mass * (z_tip - zc).powi(2);
        [self.bus_i[0] + add, self.bus_i[1] + add, self.bus_i[2]]
    }
}

fn deg(x: f64) -> f64 {
    x.to_degrees()
}

fn stability_map() {
    println!("1. Linear stability in the (k1, k3) plane, circular orbit");
    println!("   L = Lagrange (I2 > I1 > I3), D = DeBra-Delp, . = unstable");
    let n = 41;
    let (mut nl, mut nd) = (0, 0);
    for r in 0..n {
        let k3 = 1.0 - 2.0 * r as f64 / (n - 1) as f64;
        let mut row = String::new();
        for c in 0..n {
            let k1 = -1.0 + 2.0 * c as f64 / (n - 1) as f64;
            // |k| < 1 with k1 k3 < 1 is the physically realisable square
            let ch = match linear::classify(k1, k3) {
                Region::Lagrange => {
                    nl += 1;
                    'L'
                }
                Region::DeBraDelp => {
                    nd += 1;
                    'D'
                }
                Region::Unstable => '.',
            };
            row.push(ch);
        }
        if r % 2 == 0 {
            println!("   k3 {:+.2} {}", k3, row);
        }
    }
    println!("   k1 runs -1 .. +1 left to right");
    // fine-grid area fractions
    let m = 2001;
    let (mut al, mut ad) = (0usize, 0usize);
    for a in 0..m {
        for b in 0..m {
            let k1 = -1.0 + 2.0 * (a as f64 + 0.5) / m as f64;
            let k3 = -1.0 + 2.0 * (b as f64 + 0.5) / m as f64;
            match linear::classify(k1, k3) {
                Region::Lagrange => al += 1,
                Region::DeBraDelp => ad += 1,
                _ => {}
            }
        }
    }
    let tot = (m * m) as f64;
    println!(
        "   area of the k-square: Lagrange {:.4} (triangle = 0.1250), DeBra-Delp {:.4}  [coarse map cells {} / {}]\n",
        al as f64 / tot,
        ad as f64 / tot,
        nl,
        nd
    );
}

fn tongue_table() {
    println!("2. Parametric pitch tongues k2 in [lo, hi] (Floquet, |tr| > 2 inside)");
    println!("   n = 1 first-order prediction: 1/12 -+ e/8; n = 3 width should scale as e^3");
    println!("      e     n=1 lo     n=1 hi    pred lo    pred hi    n=2 width   n=3 width");
    for &e in &[0.005, 0.01, 0.02, 0.05, 0.1, 0.2] {
        let t1 = pitch::tongue(1, e).unwrap();
        let w = |n| pitch::tongue(n, e).map(|(a, b)| b - a).unwrap_or(0.0);
        // n = 2 never opens: M(k2 = 1/3) = I for every e (checked in tests)
        println!(
            "   {:5.3}  {:9.6}  {:9.6}  {:9.6}  {:9.6}  {:10.3e}  {:10.3e}",
            e,
            t1.0,
            t1.1,
            1.0 / 12.0 - e / 8.0,
            1.0 / 12.0 + e / 8.0,
            w(2),
            w(3)
        );
    }
    println!();
}

fn boom_sweep(c: &Craft, e: f64) -> f64 {
    println!("3. Boom length sweep, e = {:.4}", e);
    println!("    L[m]    I1     I2     I3     k1     k2     k3   region   wp    wr1   wr2   A_lin[deg] A_nl[deg]  |tr| per.orb");
    let mut chosen = f64::NAN;
    let mut l = 0.0;
    while l <= 6.0 + 1e-9 {
        let i = c.inertia(l);
        let b = Body { i, e };
        let (k1, k2, k3) = (b.k1(), b.k2(), b.k3());
        let reg = linear::classify(k1, k3);
        let wp = linear::pitch_freq(k2).unwrap_or(f64::NAN);
        let (wr1, wr2) = linear::roll_yaw_freqs(k1, k3).unwrap_or((f64::NAN, f64::NAN));
        let a_lin = pitch::forced_amplitude_linear(k2, e);
        let po = pitch::periodic_orbit(k2, e, [0.0, a_lin]);
        let (a_nl, tr) = match po {
            Some((x0, m)) => (pitch::peak_over_orbit(k2, e, x0), m[0] + m[3]),
            None => (f64::NAN, f64::NAN),
        };
        let reg_s = match reg {
            Region::Lagrange => "Lagr",
            Region::DeBraDelp => "DbDp",
            Region::Unstable => "unst",
        };
        println!(
            "   {:5.2}  {:5.2}  {:5.2}  {:5.2}  {:5.3}  {:5.3}  {:5.3}  {}  {:5.3} {:5.3} {:5.3}  {:9.3}  {:9.3}  {:6.3}  {}",
            l,
            i[0],
            i[1],
            i[2],
            k1,
            k2,
            k3,
            reg_s,
            wp,
            wr1,
            wr2,
            deg(a_lin.abs()),
            deg(a_nl),
            tr.abs(),
            if tr.abs() < 2.0 { "stable" } else { "UNSTABLE" }
        );
        if chosen.is_nan() && reg == Region::Lagrange && deg(a_nl) <= 2.0 && tr.abs() < 2.0 {
            chosen = l;
        }
        l += 0.25;
    }
    println!();
    chosen
}

fn refine_length(c: &Craft, e: f64, lo: f64, hi: f64, limit_deg: f64) -> f64 {
    // nonlinear periodic amplitude falls monotonically with L above k2 = 1/3
    let amp = |l: f64| {
        let b = Body { i: c.inertia(l), e };
        let k2 = b.k2();
        let (x0, _) =
            pitch::periodic_orbit(k2, e, [0.0, pitch::forced_amplitude_linear(k2, e)]).unwrap();
        deg(pitch::peak_over_orbit(k2, e, x0))
    };
    let (mut a, mut b) = (lo, hi);
    for _ in 0..60 {
        let m = 0.5 * (a + b);
        if amp(m) > limit_deg {
            a = m;
        } else {
            b = m;
        }
    }
    b
}

fn full_3d_check(c: &Craft, l: f64, e: f64) {
    let b = Body { i: c.inertia(l), e };
    let k2 = b.k2();
    let (x0, _) =
        pitch::periodic_orbit(k2, e, [0.0, pitch::forced_amplitude_linear(k2, e)]).unwrap();
    println!("5. Full nonlinear 3D check at L = {:.3} m, 50 orbits", l);
    let orbits = 50;
    let steps = 3000 * orbits;
    for &(r0, y0) in &[
        (0.0, 0.0),
        (0.5_f64.to_radians(), 0.5_f64.to_radians()),
        (3.0_f64.to_radians(), 3.0_f64.to_radians()),
    ] {
        let mut y = b.state(0.0, (r0, x0[0], y0), (0.0, x0[1] * b.nu_dot(0.0), 0.0));
        let (mut pr, mut pp, mut py) = (0.0f64, 0.0f64, 0.0f64);
        let mut pitch_err = 0.0f64;
        let f = |t: f64, s: &[f64], d: &mut [f64]| b.rhs(t, s, d);
        ode::rk4(&f, 0.0, 2.0 * PI * orbits as f64, &mut y, steps, |_, s| {
            let (r, p, yw) = b.angles(s);
            pr = pr.max(r.abs());
            pp = pp.max(p.abs());
            py = py.max(yw.abs());
        });
        // after an integer number of orbits the planar periodic pitch is back at x0
        if r0 == 0.0 {
            pitch_err = (b.angles(&y).1 - x0[0]).abs();
        }
        println!(
            "   roll0 = yaw0 = {:4.1} deg: peak roll {:6.3}, pitch {:6.3}, yaw {:6.3} deg{}",
            deg(r0),
            deg(pr),
            deg(pp),
            deg(py),
            if r0 == 0.0 {
                format!("  (pitch return error {:.1e} rad)", pitch_err)
            } else {
                String::new()
            }
        );
    }
    println!();
}

fn main() {
    let rp = RE + 600.0;
    let ra = RE + 880.0;
    let e = (ra - rp) / (ra + rp);
    let craft = Craft {
        bus_mass: 50.0,
        bus_i: [2.2, 2.5, 1.6],
        attach: 0.3,
        rod_density: 0.08,
        tip_mass: 1.5,
    };

    stability_map();
    tongue_table();
    let first_ok = boom_sweep(&craft, e);

    println!("4. Boom sizing for a 2 deg pitch libration limit");
    let l_amp = refine_length(&craft, e, first_ok - 0.25, first_ok, 2.0);
    // The amplitude spec alone lands the pitch frequency on 3/2 w0, the centre
    // of the n = 3 tongue. That tongue is only ~1e-7 wide in k2 here, but a
    // libration whose Floquet multipliers sit next to -1 has no margin for
    // inertia error or amplitude-dependent frequency shift, so also require
    // |w_p - 3/2| >= 0.05 (|tr M| <= 2 cos(0.1 pi) = 1.902).
    let k2_of = |l: f64| {
        (Body {
            i: craft.inertia(l),
            e,
        })
        .k2()
    };
    let length_for_k2 = |k2t: f64| {
        let (mut a, mut b) = (0.0, 6.0);
        for _ in 0..60 {
            let m = 0.5 * (a + b);
            if k2_of(m) < k2t {
                a = m;
            } else {
                b = m;
            }
        }
        b
    };
    let l_detune = length_for_k2(1.55f64.powi(2) / 3.0);
    let l_design = ((1.15 * l_amp).max(l_detune) * 10.0).ceil() / 10.0;
    let b = Body {
        i: craft.inertia(l_design),
        e,
    };
    let k2 = b.k2();
    let (x0, m) =
        pitch::periodic_orbit(k2, e, [0.0, pitch::forced_amplitude_linear(k2, e)]).unwrap();
    let a_nl = deg(pitch::peak_over_orbit(k2, e, x0));
    let (wr1, wr2) = linear::roll_yaw_freqs(b.k1(), b.k3()).unwrap();
    println!(
        "   amplitude limit alone:           L >= {:.3} m (w_p = {:.3} w0)",
        l_amp,
        (3.0 * k2_of(l_amp)).sqrt()
    );
    println!(
        "   1.15 x amplitude length:         L >= {:.3} m",
        1.15 * l_amp
    );
    println!("   detuning |w_p - 1.5| >= 0.05:    L >= {:.3} m", l_detune);
    println!("   design length (rounded up):      L  = {:.1} m", l_design);
    println!(
        "   k1 {:.4}  k2 {:.4}  k3 {:.4}, pitch {:.3} w0, roll/yaw {:.3} / {:.3} w0",
        b.k1(),
        k2,
        b.k3(),
        (3.0 * k2).sqrt(),
        wr1,
        wr2
    );
    println!(
        "   periodic pitch libration {:.3} deg (linear {:.3}), tr M = {:.4}",
        a_nl,
        deg(pitch::forced_amplitude_linear(k2, e).abs()),
        m[0] + m[3]
    );
    if let Some((lo, hi)) = pitch::tongue(3, e) {
        println!(
            "   n = 3 tongue at this e: k2 in [{:.7}, {:.7}] (width {:.1e}), reached at L = {:.3} m",
            lo,
            hi,
            hi - lo,
            length_for_k2(lo)
        );
    }
    // +-10 % on the bus yaw inertia, the least certain number before the boom
    // is weighed, moves k2 by:
    let mut lo_k2 = f64::MAX;
    let mut hi_k2 = f64::MIN;
    for &s in &[0.9, 1.1] {
        let mut i = craft.inertia(l_design);
        i[2] *= s;
        let kk = (i[0] - i[2]) / i[1];
        lo_k2 = lo_k2.min(kk);
        hi_k2 = hi_k2.max(kk);
    }
    println!(
        "   +-10 % bus yaw inertia: k2 {:.4} .. {:.4}, w_p {:.3} .. {:.3} w0 (stays clear of 1.5)",
        lo_k2,
        hi_k2,
        (3.0 * lo_k2).sqrt(),
        (3.0 * hi_k2).sqrt()
    );
    println!(
        "   deployment sweeps through the k2 = 1/3 forced resonance at L = {:.3} m\n",
        length_for_k2(1.0 / 3.0)
    );

    full_3d_check(&craft, l_design, e);
}
