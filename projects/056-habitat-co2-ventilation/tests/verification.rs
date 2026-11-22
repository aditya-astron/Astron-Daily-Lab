// Astron Built this AP
// Project 056 | Analog habitat systems / ECLSS | Multi-zone CO2 ventilation network for a sealed habitat
// Author: Aditya (aditya-astron)
// Date: 2026-11-22 IST
//
// Checks of the propagator, eigensolver and zone model against closed forms,
// conservation laws and an independent RK4 integration.

use habitat_co2::habitat::*;
use habitat_co2::linalg::*;

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs().max(1e-300)
}

#[test]
fn expm_matches_rotation_and_diagonal() {
    let t = 2.3;
    let mut a = Mat::zeros(2);
    a[(0, 1)] = t;
    a[(1, 0)] = -t;
    let e = expm(&a);
    assert!((e[(0, 0)] - t.cos()).abs() < 1e-14);
    assert!((e[(0, 1)] - t.sin()).abs() < 1e-14);
    assert!((e[(1, 0)] + t.sin()).abs() < 1e-14);

    let mut d = Mat::zeros(3);
    d[(0, 0)] = -40.0;
    d[(1, 1)] = 0.3;
    d[(2, 2)] = 5.0;
    let e = expm(&d);
    for (i, v) in [-40.0f64, 0.3, 5.0].iter().enumerate() {
        assert!(rel(e[(i, i)], v.exp()) < 1e-13);
    }
}

#[test]
fn jacobi_matches_tridiagonal_spectrum() {
    // tridiag(-1, 2, -1) of order n has eigenvalues 2 - 2 cos(k pi / (n+1)).
    let n = 6;
    let mut a = Mat::zeros(n);
    for i in 0..n {
        a[(i, i)] = 2.0;
        if i + 1 < n {
            a[(i, i + 1)] = -1.0;
            a[(i + 1, i)] = -1.0;
        }
    }
    let ev = sym_eigenvalues(&a);
    for k in 1..=n {
        let exact = 2.0 - 2.0 * (k as f64 * std::f64::consts::PI / (n as f64 + 1.0)).cos();
        assert!((ev[k - 1] - exact).abs() < 1e-13);
    }
}

#[test]
fn unit_conversion_hand_check() {
    // 1 kg/m^3 of CO2 at 295.15 K: p = rho R T / M = 55 760 Pa = 418.2 mmHg.
    assert!((mmhg(1.0) - 418.24).abs() < 0.01);
    assert!(rel(density(3.0), 3.0 / mmhg(1.0)) < 1e-15);
    // 1.04 kg/d over 16 h awake + 8 h at 0.7 rate.
    assert!(rel(16.0 * awake_rate() + 8.0 * sleep_rate(), CO2_PER_DAY) < 1e-15);
}

#[test]
fn single_zone_exponential_approach() {
    // One scrubbed zone: c(t) = c_inf + (c0 - c_inf) exp(-eta Qs t / V).
    let h = Habitat { names: vec!["cabin"], vol: vec![50.0], links: vec![], scrub_zone: 0, scrub_flow: 40.0, scrub_eff: 0.9 };
    let s = [0.1];
    let c_inf = 0.1 / 36.0;
    let c0 = [0.001];
    for t in [0.1, 1.0, 5.0] {
        let exact = c_inf + (c0[0] - c_inf) * (-36.0 * t / 50.0f64).exp();
        assert!(rel(h.advance(&c0, &s, t)[0], exact) < 1e-13);
    }
    assert!(rel(h.steady(&s)[0], c_inf) < 1e-14);
}

#[test]
fn two_zone_eigenvalues_from_the_quadratic() {
    // Zones a (scrubbed, removal R) and b with exchange Q:
    // det [[-(Q+R)/Va - l, Q/Va], [Q/Vb, -Q/Vb - l]] = 0.
    let (va, vb, q, r) = (80.0, 30.0, 25.0, 40.0);
    let h = Habitat { names: vec!["a", "b"], vol: vec![va, vb], links: vec![(0, 1, q)], scrub_zone: 0, scrub_flow: r, scrub_eff: 1.0 };
    let tr = (q + r) / va + q / vb;
    let det = r * q / (va * vb);
    let disc = (tr * tr - 4.0 * det).sqrt();
    let exact = [(tr - disc) / 2.0, (tr + disc) / 2.0];
    let rates = h.decay_rates();
    assert!(rel(rates[0], exact[0]) < 1e-13);
    assert!(rel(rates[1], exact[1]) < 1e-13);
}

#[test]
fn hub_closed_form_matches_linear_solve() {
    for q in [10.0, 38.0, 60.0, 200.0] {
        let h = reference_habitat(q);
        for src in [day_source(), sleep_source()] {
            let a = h.steady(&src);
            let b = hub_steady_closed_form(&h, NODE, &src);
            for k in 0..h.n() {
                assert!(rel(a[k], b[k]) < 1e-12, "q={q} zone {k}: {} vs {}", a[k], b[k]);
            }
        }
    }
}

#[test]
fn spectrum_invariants_and_rayleigh_bound() {
    let h = reference_habitat(60.0);
    let rates = h.decay_rates();
    let a = h.system();
    // Sum of eigenvalues equals the trace of A = V^-1 L.
    assert!(rel(rates.iter().sum::<f64>(), -a.trace()) < 1e-13);
    // Rayleigh quotient of V^-1/2 L V^-1/2 at V^1/2 * 1 gives eta Qs / V_total,
    // so the slowest decay rate can be no faster than the well-mixed one.
    let vt: f64 = h.vol.iter().sum();
    assert!(rates[0] <= h.scrub_eff * h.scrub_flow / vt);
    assert!(rates.iter().all(|&r| r > 0.0));
}

#[test]
fn exact_map_matches_rk4() {
    let h = reference_habitat(60.0);
    let c0 = vec![0.003, 0.004, 0.0035, 0.002, 0.006];
    for (src, t) in [(day_source(), 16.0), (sleep_source(), 8.0)] {
        let exact = h.advance(&c0, &src, t);
        // Fastest mode is 0.07 h; dt = 2e-3 h keeps RK4 error near 1e-13.
        let steps = (t / 2e-3) as usize;
        let num = h.rk4(&c0, &src, t / steps as f64, steps);
        for k in 0..h.n() {
            assert!(rel(exact[k], num[k]) < 1e-10, "zone {k}: {} vs {}", exact[k], num[k]);
        }
    }
}

#[test]
fn closed_habitat_conserves_mass() {
    // No scrubber: total CO2 mass sum V c rises by exactly S t.
    let mut h = reference_habitat(60.0);
    h.scrub_eff = 0.0;
    let c0 = vec![0.003; 5];
    let src = day_source();
    let t = 5.0;
    let c = h.advance(&c0, &src, t);
    let m0: f64 = h.vol.iter().zip(&c0).map(|(v, c)| v * c).sum();
    let m1: f64 = h.vol.iter().zip(&c).map(|(v, c)| v * c).sum();
    assert!(rel(m1 - m0, src.iter().sum::<f64>() * t) < 1e-12);
}

#[test]
fn periodic_state_is_periodic_and_attracting() {
    let h = reference_habitat(60.0);
    let sched = reference_schedule();
    let x0 = h.periodic_start(&sched);
    let mut x = x0.clone();
    for ph in &sched {
        x = h.advance(&x, &ph.source, ph.hours);
    }
    for k in 0..h.n() {
        assert!(rel(x[k], x0[k]) < 1e-12);
    }
    // Brute force: 30 days from a clean habitat, RK4 at dt = 0.004 h.
    // Slowest time constant is 7.8 h, so 720 h leaves e^-92 of the start.
    let mut y = vec![0.0; h.n()];
    for _ in 0..30 {
        for ph in &sched {
            let steps = (ph.hours / 0.004).round() as usize;
            y = h.rk4(&y, &ph.source, ph.hours / steps as f64, steps);
        }
    }
    for k in 0..h.n() {
        assert!(rel(y[k], x0[k]) < 1e-9, "zone {k}: {} vs {}", y[k], x0[k]);
    }
}

#[test]
fn periodic_mean_obeys_scrubber_mass_balance() {
    // Over one period the hull holds the same CO2 at both ends, so the
    // scrubber removes the daily production: eta Qs <c_Hab> = <S>.
    let h = reference_habitat(45.0);
    let sched = reference_schedule();
    let mut x = h.periodic_start(&sched);
    let mut integral = 0.0;
    for ph in &sched {
        let steps = (ph.hours * 120.0) as usize;
        let dt = ph.hours / steps as f64;
        let (phi, g) = h.flow_map(&ph.source, dt);
        for _ in 0..steps {
            let nx: Vec<f64> = phi.mul_vec(&x).iter().zip(&g).map(|(u, v)| u + v).collect();
            integral += 0.5 * dt * (x[HAB] + nx[HAB]);
            x = nx;
        }
    }
    let removed = h.scrub_eff * h.scrub_flow * integral;
    let made: f64 = sched.iter().map(|p| p.hours * p.source.iter().sum::<f64>()).sum();
    assert!(rel(removed, made) < 1e-6);
}

#[test]
fn time_to_limit_brackets_crossing() {
    let h = reference_habitat(3.0);
    let c0 = vec![0.0028, 0.0037, 0.0037, 0.0037, 0.0039];
    let lim = density(3.0);
    let t = h.time_to_limit(&c0, &sleep_source(), QUARTERS, lim, 24.0).unwrap();
    assert!(rel(h.advance(&c0, &sleep_source(), t)[QUARTERS], lim) < 1e-12);
    assert!(h.advance(&c0, &sleep_source(), t - 1e-3)[QUARTERS] < lim);
    // Sealed-room estimate (no exchange): t ~ V dc / S, slightly shorter.
    let sealed = h.vol[QUARTERS] * (lim - c0[QUARTERS]) / sleep_source()[QUARTERS];
    assert!(t > sealed && t < 1.1 * sealed);
}

#[test]
fn peak_search_finds_interior_maximum() {
    let h = reference_habitat(60.0);
    let sched = reference_schedule();
    let x0 = h.periodic_start(&sched);
    let xs = h.advance(&x0, &sched[0].source, 16.0);
    let (t, p) = h.phase_peak(&xs, &sched[1], QUARTERS);
    assert!(t > 0.0 && t < 8.0);
    for dt in [-1e-3, 1e-3] {
        assert!(h.advance(&xs, &sched[1].source, t + dt)[QUARTERS] <= p);
    }
}
