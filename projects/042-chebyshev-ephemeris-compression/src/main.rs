// Astron Built this AP
// Project 042 | Numerical methods / Navigation | Chebyshev ephemeris compression for an eccentric orbit
// Author: Aditya (aditya-astron)
// Date: 2026-11-16 IST

use cheb_ephemeris::*;

fn molniya() -> Orbit {
    Orbit {
        a: 26_554.0,
        e: 0.74,
        inc: 63.4f64.to_radians(),
        raan: 40.0f64.to_radians(),
        argp: 270.0f64.to_radians(),
        mu: MU_EARTH,
    }
}

fn main() {
    let orb = molniya();
    let n_orb = orb.mean_motion();
    println!("Molniya orbit: a = {:.0} km, e = {:.2}, period = {:.1} min", orb.a, orb.e, orb.period() / 60.0);
    println!("Kepler branch points at M* = 2 pi m +/- i eta, eta = {:.4} rad", orb.eta());

    // 1. Convergence rate on a single segment vs the Bernstein prediction.
    println!("\n1) Coefficient decay of x(t) on one segment, measured vs predicted rho");
    println!("   centre      half-width   rho pred   rho meas");
    let pi = std::f64::consts::PI;
    for (label, mc, hms) in [("perigee", 0.0, [0.15, 0.3, 0.6]), ("M = 90 deg", 0.5 * pi, [0.6, 1.0, 1.5]), ("apogee", pi, [1.0, 1.5, 2.5])] {
        for hm in hms {
            let (t0, t1) = ((mc - hm) / n_orb, (mc + hm) / n_orb);
            let s = Segment::fit(&orb, t0, t1, 80);
            let pred = segment_rho(&orb, t0, t1);
            let c0 = s.coef[0][0].abs().max(s.coef[0][1].abs());
            // Fit the envelope over the range clearly above round-off.
            let kmax = (1..s.coef[0].len()).filter(|&k| s.coef[0][k].abs() > 1e-11 * c0).max().unwrap_or(10);
            let k_lo = 6.max(kmax / 3);
            let meas = if kmax > k_lo + 6 { measured_rho(&s.coef[0], k_lo, kmax - 3) } else { f64::NAN };
            println!("   {:<11} {:>6.2} rad    {:>7.3}    {:>7.3}  (k {}..{})", label, hm, pred, meas, k_lo, kmax - 3);
        }
    }

    // 2. Layout trade at 1 mm. Doubles counted are stored coefficients only.
    let tol = 1e-6;
    println!("\n2) Doubles per orbit for 1 mm max position error (checked on {} points per segment)", CHECK_SAMPLES);
    println!("    N | uniform seg/dbl | bisection seg/dbl | equal-rho seg/dbl   rho*  | max vel err (equal-rho)");
    for n in [7usize, 10, 13, 18, 24] {
        let u = uniform_layout(&orb, n, tol, 4000).unwrap();
        let a = adaptive_layout(&orb, n, tol);
        let (rs, q) = tuned_equal_rho_layout(&orb, n, tol);
        let ev = q.iter().map(|s| s.max_error(&orb, CHECK_SAMPLES).1).fold(0.0, f64::max);
        let d = |k: usize| k * 3 * (n + 1);
        println!(
            "   {:>2} | {:>5} / {:>5}    | {:>5} / {:>5}      | {:>5} / {:>5}    {:>6.2} | {:.2} mm/s",
            n, u.len(), d(u.len()), a.len(), d(a.len()), q.len(), d(q.len()), rs, ev * 1e6
        );
    }
    let kh = hermite_knots_for(&orb, tol);
    println!("   cubic Hermite on (r, v) at uniform knots: {} knots, {} doubles", kh, kh * 6);

    // 3. The N = 13 equal-rho layout, segment by segment.
    let (_, q) = tuned_equal_rho_layout(&orb, 13, tol);
    println!("\n3) N = 13 equal-rho layout (times from perigee)");
    println!("   seg    t0 [min]   span [min]   rho    max |dr| [mm]");
    for (i, s) in q.iter().enumerate() {
        let (ep, _) = s.max_error(&orb, 2000);
        println!(
            "   {:>3}  {:>9.2}   {:>9.2}   {:>6.2}   {:>8.3}",
            i, s.t0 / 60.0, (s.t1 - s.t0) / 60.0, segment_rho(&orb, s.t0, s.t1), ep * 1e6
        );
    }
}
