// Astron Built this AP
// Project 016 | Orbital mechanics | Sun-synchronous repeat ground track design
// Author: Aditya (aditya-astron)
// Date: 2026-11-05 IST

use sso_repeat_track::*;
use std::f64::consts::PI;

fn close(x: f64, y: f64, tol: f64) -> bool {
    (x - y).abs() <= tol
}

#[test]
fn design_is_sun_synchronous_and_closes_exactly() {
    for (n, d) in [(233, 16), (175, 12), (143, 10), (29, 2)] {
        let o = solve_repeat_sso(n, d).unwrap();
        let r = secular_rates(o.a, o.inc);
        assert!(close(r.raan_dot / SUN_RATE, 1.0, 1e-12), "{n}/{d} not SSO");
        // N nodal periods must span exactly D nodal days.
        let lhs = n as f64 * o.nodal_period;
        let rhs = d as f64 * o.nodal_day;
        assert!(close(lhs / rhs, 1.0, 1e-11), "{n}/{d}: {lhs} vs {rhs}");
    }
}

#[test]
fn sso_nodal_day_is_the_mean_solar_day() {
    // omega_E - Sun rate = 2*pi / 86400 s up to the sidereal/tropical constants.
    let nd = nodal_day(SUN_RATE);
    assert!(close(nd, 86_400.0, 0.05), "nodal day {nd}");
}

#[test]
fn matches_published_mission_orbits() {
    // Nominal altitude and inclination from mission documentation. Published
    // altitudes are not all defined as a - Re and real orbits are tuned with
    // higher-order terms, so a few km and ~0.1 deg is the expected agreement.
    let cases = [
        ("Landsat-8", 233, 16, 705.0, 98.2),
        ("Sentinel-1", 175, 12, 693.0, 98.18),
        ("Sentinel-2", 143, 10, 786.0, 98.62),
    ];
    for (name, n, d, h, i) in cases {
        let o = solve_repeat_sso(n, d).unwrap();
        assert!(close(o.altitude(), h, 6.0), "{name}: h = {}", o.altitude());
        assert!(close(o.inc.to_degrees(), i, 0.1), "{name}: i = {}", o.inc.to_degrees());
    }
}

#[test]
fn j2_correction_has_the_expected_size() {
    // Hand estimate: without J2, n_kep = (N/D) * 2*pi/86400 gives a Keplerian a.
    // With J2 and i ~ 98 deg, dM/dt + dw/dt = n [1 + 3/4 k (8 cos^2 i - 2)] is
    // slower than n by about 1.5 k, so the J2 orbit sits lower by about
    // a * (2/3) * 1.5 k = a k (~ 8.8 km for Landsat).
    let o = solve_repeat_sso(233, 16).unwrap();
    let n_kep = 233.0 / 16.0 * 2.0 * PI / o.nodal_day;
    let a_kep = (MU / (n_kep * n_kep)).cbrt();
    let k = J2 * (RE / o.a).powi(2);
    let c2 = o.inc.cos().powi(2);
    let expected = o.a * (2.0 / 3.0) * 0.75 * k * (2.0 - 8.0 * c2);
    assert!(a_kep > o.a);
    assert!(close(a_kep - o.a, expected, 0.05), "{} vs {}", a_kep - o.a, expected);
}

/// Brute-force sub-cycle: walk revolution by revolution and find the first
/// rev whose node lands one track spacing east or west of the reference.
fn subcycle_brute(n: u32, d: u32) -> (u32, bool) {
    for k in 1..=n {
        let r = (k * d) % n; // node is r spacings west of the reference
        if r == 1 || r == n - 1 {
            let days = (k as f64 * d as f64 / n as f64).round() as u32;
            return (days, r == n - 1);
        }
    }
    unreachable!()
}

#[test]
fn subcycle_matches_brute_force_and_landsat() {
    // Landsat WRS-2: adjacent paths are imaged 7 days apart.
    assert_eq!(subcycle(233, 16).0, 7);
    for d in 2..=30u32 {
        for n in (13 * d)..=(16 * d) {
            let mut g = (n, d);
            while g.1 != 0 {
                g = (g.1, g.0 % g.1);
            }
            if g.0 != 1 {
                continue;
            }
            assert_eq!(subcycle(n, d), subcycle_brute(n, d), "{n}/{d}");
        }
    }
}

#[test]
fn enumeration_is_complete_and_valid() {
    let list = enumerate_repeat_orbits(600.0, 800.0, 7);
    for o in &list {
        assert!(o.altitude() >= 600.0 && o.altitude() <= 800.0);
        let lhs = o.revs as f64 * o.nodal_period;
        assert!(close(lhs / (o.days as f64 * o.nodal_day), 1.0, 1e-11));
    }
    // Known members of the band.
    for (n, d) in [(29, 2), (44, 3), (73, 5), (102, 7)] {
        assert!(list.iter().any(|o| o.revs == n && o.days == d), "{n}/{d} missing");
    }
    // Daily repeat 14 or 15 rev/day sits outside 600-800 km.
    assert!(!list.iter().any(|o| o.days == 1));
}

#[test]
fn short_period_offset_at_node_matches_first_order_theory() {
    // Osculating minus mean a for a circular orbit: 3/2 J2 Re^2/a sin^2 i cos 2u,
    // evaluated at the ascending node (u = 0).
    let o = solve_repeat_sso(233, 16).unwrap();
    let (a_osc, _) = mean_to_osculating_at_node(o.a, o.inc, o.nodal_period);
    let expected = 1.5 * J2 * RE * RE / o.a * o.inc.sin().powi(2);
    assert!(close(a_osc - o.a, expected, 0.1), "{} vs {}", a_osc - o.a, expected);
}

#[test]
fn numerical_propagation_reproduces_the_design() {
    let o = solve_repeat_sso(233, 16).unwrap();
    let chk = verify_by_propagation(&o, 5.0);
    // First-order theory vs full J2 dynamics: period to ~1e-6, node rate to
    // within the O(J2^2) terms (~0.1%), and the track must close to well under
    // 1% of the 172 km track spacing after 16 days.
    assert!(close(chk.nodal_period / o.nodal_period, 1.0, 1e-5));
    assert!(close(chk.raan_rate / SUN_RATE, 1.0, 2e-3));
    assert!(chk.closure_km.abs() < 1.72, "closure {} km", chk.closure_km);
    // Step-size convergence: halving the RK4 step changes nothing material.
    let fine = verify_by_propagation(&o, 2.5);
    assert!(close(fine.closure_km, chk.closure_km, 0.01));
}
