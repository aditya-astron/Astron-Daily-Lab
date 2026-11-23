// Astron Built this AP
// Project 061 | Radiation shielding | Sector-based TID analysis and spot-shield sizing for an avionics box
// Author: Aditya (aditya-astron)
// Date: 2026-11-23 IST
//
// Verification of the special functions, ray-body chords, the sector
// integrator against closed forms, and the sizing searches.

use sector_tid::dose::*;
use sector_tid::geometry::*;
use sector_tid::model::*;
use sector_tid::sector::*;

fn rel(a: f64, b: f64) -> f64 {
    (a / b - 1.0).abs()
}

#[test]
fn exponential_integrals_match_tabulated_values() {
    // Abramowitz & Stegun table 5.1 / DLMF values.
    let e1_ref = [
        (0.1, 1.822_923_958_419_390_7),
        (0.5, 0.559_773_594_776_160_8),
        (1.0, 0.219_383_934_395_520_3),
        (2.0, 0.048_900_510_708_061_1),
        (10.0, 4.156_968_929_685_324e-6),
    ];
    for (x, v) in e1_ref {
        assert!(rel(e1(x), v) < 1e-13, "E1({x}) = {} vs {v}", e1(x));
    }
    assert!(rel(e2(1.0), 0.148_495_506_775_922_0) < 1e-13);
}

#[test]
fn e2_matches_its_defining_integral_across_series_and_fraction_branches() {
    for x in [0.02, 0.3, 0.999, 1.001, 2.5, 7.0, 20.0] {
        let q = adaptive_simpson(&|mu: f64| if mu == 0.0 { 0.0 } else { (-x / mu).exp() }, 0.0, 1.0, 1e-14 * (-x).exp());
        assert!(rel(e2(x), q) < 1e-10, "x = {x}: {} vs {q}", e2(x));
    }
}

#[test]
fn box_and_sphere_chords() {
    let lo = [0.0, 0.0, 0.0];
    let hi = [2.0, 3.0, 4.0];
    // From inside along +x: distance to the far face.
    assert!((chord_box(&[0.5, 1.0, 1.0], &[1.0, 0.0, 0.0], &lo, &hi) - 1.5).abs() < 1e-15);
    // From outside, straight through.
    assert!((chord_box(&[-5.0, 1.0, 1.0], &[1.0, 0.0, 0.0], &lo, &hi) - 2.0).abs() < 1e-15);
    // Pointing away, and parallel outside the slab.
    assert_eq!(chord_box(&[-5.0, 1.0, 1.0], &[-1.0, 0.0, 0.0], &lo, &hi), 0.0);
    assert_eq!(chord_box(&[-5.0, 5.0, 1.0], &[1.0, 0.0, 0.0], &lo, &hi), 0.0);
    // Space diagonal of a unit cube.
    let s = 1.0 / 3.0_f64.sqrt();
    assert!((chord_box(&[-1.0, -1.0, -1.0], &[s, s, s], &[0.0; 3], &[1.0; 3]) - 3.0_f64.sqrt()).abs() < 1e-14);

    assert!((chord_sphere(&[0.0; 3], &[0.0, 0.0, 1.0], &[0.0; 3], 2.0) - 2.0).abs() < 1e-15);
    assert!((chord_sphere(&[0.0, 0.0, -5.0], &[0.0, 0.0, 1.0], &[0.0; 3], 2.0) - 4.0).abs() < 1e-14);
    // Chord at impact parameter 1 through radius 2: 2 sqrt(3).
    assert!((chord_sphere(&[1.0, 0.0, -5.0], &[0.0, 0.0, 1.0], &[0.0; 3], 2.0) - 2.0 * 3.0_f64.sqrt()).abs() < 1e-14);
    assert_eq!(chord_sphere(&[3.0, 0.0, -5.0], &[0.0, 0.0, 1.0], &[0.0; 3], 2.0), 0.0);
}

#[test]
fn hollow_bodies_subtract_cavity_exactly() {
    let mut s = Scene::new();
    s.add_box_shell([-10.0; 3], [10.0; 3], 0.5, 2.0);
    let t = s.areal_density(&[0.0; 3], &[0.0, 1.0, 0.0]);
    assert!((t - 1.0).abs() < 1e-13);
    let t = s.areal_density(&[0.0, -20.0, 0.0], &[0.0, 1.0, 0.0]);
    assert!((t - 2.0).abs() < 1e-13);
    let m = 2.0 * (20.0_f64.powi(3) - 19.0_f64.powi(3));
    assert!(rel(s.mass(), m) < 1e-14);

    let mut s = Scene::new();
    s.add_sphere_shell([0.0; 3], 5.0, 1.0, 3.0);
    assert!(rel(s.mass(), 3.0 * 4.0 / 3.0 * std::f64::consts::PI * (125.0 - 64.0)) < 1e-14);
}

#[test]
fn fibonacci_lattice_is_unit_and_balanced() {
    let n = 20_000;
    let d = fibonacci_dirs(n);
    let mut m = [0.0; 3];
    let mut zz = 0.0;
    for u in &d {
        assert!(((u[0] * u[0] + u[1] * u[1] + u[2] * u[2]) - 1.0).abs() < 1e-14);
        for k in 0..3 {
            m[k] += u[k] / n as f64;
        }
        zz += u[2] * u[2] / n as f64;
    }
    assert!(m.iter().all(|x| x.abs() < 1e-3), "{m:?}");
    // Midpoint rule in z: error in <z^2> is exactly 1/(3 N^2).
    assert!((zz - 1.0 / 3.0 + 1.0 / (3.0 * (n * n) as f64)).abs() < 1e-14);
}

#[test]
fn sector_dose_matches_closed_forms() {
    let dd = DoseDepth::geo_15yr_illustrative();
    let dirs = fibonacci_dirs(200_000);

    let mut s = Scene::new();
    s.add_sphere_shell([0.0; 3], 20.0, 0.3, AL_RHO);
    assert!(rel(sector_dose(&s, &[0.0; 3], &dd, &dirs), dd.eval(0.3 * AL_RHO)) < 1e-12);
    let q = offset_shell_reference(&dd, 19.7, 20.0, AL_RHO, 12.0, 1e-12);
    assert!(rel(sector_dose(&s, &[0.0, 0.0, 12.0], &dd, &dirs), q) < 1e-9);

    for w in [0.1, 0.54, 2.0] {
        let mut s = Scene::new();
        s.add_box([-1e6, -1e6, 1.0], [1e6, 1e6, 1.0 + w], 1.0);
        s.add_box([-1e6, -1e6, -1.0 - w], [1e6, 1e6, -1.0], 1.0);
        let got = sector_dose(&s, &[0.0; 3], &dd, &dirs);
        assert!(rel(got, dd.two_sided_slab(w)) < 1e-9, "w = {w}");
    }

    let mut s = Scene::new();
    s.add_box_shell([-10.0; 3], [10.0; 3], 0.2, AL_RHO);
    let q = cube_centre_reference(&dd, 0.2 * AL_RHO, 1e-13);
    assert!(rel(sector_dose(&s, &[0.0; 3], &dd, &dirs), q) < 1e-6);
}

#[test]
fn cube_quadrature_covers_the_full_sphere() {
    // A flat dose-depth curve must integrate to itself.
    let flat = DoseDepth { terms: vec![(1.0, 1e30)] };
    assert!((cube_centre_reference(&flat, 1.0, 1e-14) - 1.0).abs() < 1e-12);
}

#[test]
fn equivalent_thickness_inverts_the_curve() {
    let dd = DoseDepth::geo_15yr_illustrative();
    for t in [0.05, 0.5, 1.6, 4.0] {
        assert!(rel(dd.equivalent_thickness(dd.eval(t)), t) < 1e-10);
    }
    assert!((mm_al_to_areal(1.0) - 0.27).abs() < 1e-15);
}

#[test]
fn fibonacci_beats_random_sampling_in_the_bus_model() {
    let dd = DoseDepth::geo_15yr_illustrative();
    let mut s = bus_scene();
    s.extend(&avionics_box(0.0, 0.0));
    let p = parts()[1].point();
    let r = sector_dose(&s, &p, &dd, &fibonacci_dirs(2_000_000));
    let e10k = rel(sector_dose(&s, &p, &dd, &fibonacci_dirs(10_000)), r);
    let e100k = rel(sector_dose(&s, &p, &dd, &fibonacci_dirs(100_000)), r);
    let er = rel(sector_dose(&s, &p, &dd, &random_dirs(100_000, 7)), r);
    assert!(e100k < e10k && e100k < 1e-4 && e100k < er, "{e10k} {e100k} {er}");
}

#[test]
fn spot_can_geometry_and_mass() {
    let p = &parts()[1];
    let w = 0.07;
    let can = spot_can(p, w, false);
    let (c, h) = (CAN_HALF, CAN_HEIGHT);
    let m = AL_RHO * ((2.0 * c + 2.0 * w).powi(2) * w + 2.0 * h * w * (2.0 * c + 2.0 * w) + 2.0 * h * w * 2.0 * c);
    assert!(rel(can.mass(), m) < 1e-13);
    // Looking straight up through the lid adds exactly w of aluminium;
    // looking sideways through a wall adds w as well.
    let pt = p.point();
    assert!((can.areal_density(&pt, &[1.0, 0.0, 0.0]) - w * AL_RHO).abs() < 1e-13);
    assert!((can.areal_density(&pt, &[0.0, 1.0, 0.0]) - w * AL_RHO).abs() < 1e-13);
    // Straight down sees nothing without the back plate.
    assert_eq!(can.areal_density(&pt, &[-1.0, 0.0, 0.0]), 0.0);
    let back = spot_can(p, w, true);
    assert!((back.areal_density(&pt, &[-1.0, 0.0, 0.0]) - w * AL_RHO).abs() < 1e-13);
    // Cans and doublers stay clear of cards and walls: no overlap means
    // total mass is additive.
    let mut s = avionics_box(0.0, 0.5);
    let m0 = s.mass();
    s.extend(&can);
    assert!(rel(s.mass(), m0 + m) < 1e-13);
}

#[test]
fn exit_face_along_axes() {
    let p = parts()[1].point();
    let axes = [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, -1.0]];
    for (f, u) in axes.iter().enumerate() {
        assert_eq!(exit_face(&p, u), f);
    }
}

#[test]
fn spot_can_sizing_lands_on_the_requirement() {
    let dd = DoseDepth::geo_15yr_illustrative();
    let dirs = fibonacci_dirs(50_000);
    let p = &parts()[1];
    let pt = p.point();
    let dose = |w: f64| {
        let mut s = bus_scene();
        s.extend(&avionics_box(0.0, 0.0));
        s.extend(&spot_can(p, w, false));
        sector_dose(&s, &pt, &dd, &dirs)
    };
    let mut prev = f64::INFINITY;
    for w in [0.0, 0.02, 0.05, 0.1, 0.2] {
        let d = dose(w);
        assert!(d < prev);
        prev = d;
    }
    let target = p.capability / 2.0;
    let w = bisect_decreasing(dose, target, 0.0, 1.0, 1e-6).unwrap();
    assert!(rel(dose(w), target) < 1e-4);
}

#[test]
fn closed_six_sided_can_is_bounded_by_the_uniform_shell_estimate() {
    // A closed can adds at least w along every ray (w / |cos| through a
    // face), so it must do better than adding exactly w to every ray, and
    // the open five-sided can must do worse than the closed one.
    let dd = DoseDepth::geo_15yr_illustrative();
    let dirs = fibonacci_dirs(50_000);
    let p = &parts()[1];
    let pt = p.point();
    let w = 0.06;
    let mut s0 = bus_scene();
    s0.extend(&avionics_box(0.0, 0.0));
    let base = thickness_map(&s0, &pt, &dirs);
    let shifted: Vec<f64> = base.iter().map(|t| t + w * AL_RHO).collect();
    let with = |back: bool| {
        let mut s = s0.clone();
        s.extend(&spot_can(p, w, back));
        sector_dose(&s, &pt, &dd, &dirs)
    };
    let (open, closed) = (with(false), with(true));
    assert!(closed < fold(&dd, &shifted));
    assert!(closed < open && open < fold(&dd, &base));
}
