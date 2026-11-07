// Astron Built this AP
// Project 020 | Spacecraft power systems | LEO Li-ion battery sizing with a 2RC equivalent-circuit model
// Author: Aditya (aditya-astron)
// Date: 2026-11-07 IST

use leo_battery::*;

fn sys(np: usize) -> PowerSystem {
    PowerSystem {
        cell: Cell::nmc_18650(),
        ns: 8,
        np,
        p_load: 180.0,
        p_array: 400.0,
        eta_bcr: 0.95,
        eta_bdr: 0.95,
        v_cell_max: 4.10,
        i_chg_max: 1.6,
        alt_m: 500e3,
        beta_rad: 0.0,
        dt: 1.0,
    }
}

#[test]
fn orbit_period_and_eclipse_match_hand_values() {
    // a = 6878.137 km: T = 2 pi sqrt(a^3/mu) = 5677.0 s.
    // cos(pi f) = sqrt(500^2 + 2*6378.137*500)/6878.137 = 0.374304 -> f = 0.377882 (35.75 min).
    let t = orbit_period(500e3);
    assert!((t - 5677.0).abs() < 0.5, "T = {t}");
    let f = eclipse_fraction(500e3, 0.0);
    let f_ref = (6_628_137.0f64.sqrt() / 6878.137).acos() / std::f64::consts::PI;
    assert!((f - f_ref).abs() < 1e-9);
    assert!((f - 0.377882).abs() < 1e-6, "f = {f}");
}

#[test]
fn eclipse_vanishes_at_critical_beta() {
    let b = beta_no_eclipse(500e3);
    assert!((b.to_degrees() - 68.02).abs() < 0.01);
    assert!(eclipse_fraction(500e3, b - 1e-3) > 0.0);
    assert_eq!(eclipse_fraction(500e3, b + 1e-3), 0.0);
    // Eclipse shrinks monotonically with beta.
    let mut last = 1.0;
    for k in 0..68 {
        let f = eclipse_fraction(500e3, (k as f64).to_radians());
        assert!(f <= last);
        last = f;
    }
}

#[test]
fn rc_branches_match_analytic_step_response_for_any_step() {
    // Under constant current Vk(t) = I Rk (1 - exp(-t/tau_k)); the ZOH update
    // is exact, so the result must not depend on the step size.
    let c = Cell::nmc_18650();
    let i = 3.0;
    for dt in [0.1, 1.0, 7.0, 45.0] {
        let mut s = CellState::rested(0.8);
        let n = (900.0 / dt) as usize;
        for _ in 0..n {
            s = step_cell(&c, &s, i, dt);
        }
        let t = n as f64 * dt;
        let v1 = i * c.r1 * (1.0 - (-t / (c.r1 * c.c1)).exp());
        let v2 = i * c.r2 * (1.0 - (-t / (c.r2 * c.c2)).exp());
        assert!((s.v1 - v1).abs() < 1e-12, "dt {dt}: v1 {} vs {}", s.v1, v1);
        assert!((s.v2 - v2).abs() < 1e-12, "dt {dt}: v2 {} vs {}", s.v2, v2);
        // Coulomb counting.
        assert!((s.soc - (0.8 - i * t / c.q_as())).abs() < 1e-12);
    }
    // Relaxation at zero current decays with tau.
    let s = CellState { soc: 0.5, v1: 0.04, v2: 0.06 };
    let r = step_cell(&c, &s, 0.0, 600.0);
    assert!((r.v2 - 0.06 * (-1.0f64).exp()).abs() < 1e-12);
    assert!((r.v1 - 0.04 * (-30.0f64).exp()).abs() < 1e-15);
}

#[test]
fn ocv_integral_matches_numerical_quadrature() {
    let c = Cell::nmc_18650();
    for &s in &[0.03, 0.27, 0.5, 0.91, 1.0] {
        let n = 200_000;
        let h = s / n as f64;
        let mut w = 0.0;
        for k in 0..n {
            w += c.ocv((k as f64 + 0.5) * h) * h;
        }
        assert!((c.ocv_integral(s) - w).abs() < 1e-9, "soc {s}");
    }
    // 4.10 V lies 20% of the way from 4.08 V (0.9) to 4.18 V (1.0).
    assert!((c.soc_at_ocv(4.10) - 0.92).abs() < 1e-12);
    assert!((c.ocv(c.soc_at_ocv(3.80)) - 3.80).abs() < 1e-12);
}

#[test]
fn constant_power_currents_deliver_the_requested_power() {
    let c = Cell::nmc_18650();
    let s = CellState { soc: 0.7, v1: 0.01, v2: -0.02 };
    for p in [0.5, 5.0, 25.0, 60.0] {
        let i = current_for_discharge_power(&c, &s, p).unwrap();
        assert!((terminal_voltage(&c, &s, i) * i - p).abs() < 1e-10);
        let ic = current_for_charge_power(&c, &s, p);
        assert!((terminal_voltage(&c, &s, -ic) * ic - p).abs() < 1e-10);
    }
    // Beyond E^2 / (4 R0) no operating point exists.
    let e = c.ocv(0.7) - 0.01 + 0.02;
    assert!(current_for_discharge_power(&c, &s, e * e / (4.0 * c.r0) * 1.001).is_none());
}

#[test]
fn ideal_cell_eclipse_discharge_matches_energy_integral() {
    // With no resistance, constant power P drains chemical energy at rate P:
    //     Q (W(s0) - W(s1)) = P t_e,
    // which fixes the end-of-eclipse SOC exactly.
    let mut s = sys(5);
    s.cell.r0 = 0.0;
    s.cell.r1 = 1e-9;
    s.cell.c1 = 1.0;
    s.cell.r2 = 1e-9;
    s.cell.c2 = 1.0;
    s.p_array = s.p_load; // no surplus in sunlight, battery idles
    s.dt = 0.25;
    let s0 = 0.9;
    let (_, st) = s.run_orbit(CellState::rested(s0));
    let p_cell = s.p_load / s.eta_bdr / s.n_cells();
    let w1 = s.cell.ocv_integral(s0) - p_cell * s.t_eclipse() / s.cell.q_as();
    let (mut lo, mut hi) = (0.0, s0);
    for _ in 0..80 {
        let m = 0.5 * (lo + hi);
        if s.cell.ocv_integral(m) < w1 { lo = m } else { hi = m }
    }
    assert!((st.soc_min - lo).abs() < 2e-5, "sim {} vs exact {}", st.soc_min, lo);
    assert!((st.soc_end - st.soc_min).abs() < 1e-12);
    assert!(st.e_loss_wh < 1e-6);
}

#[test]
fn steady_state_is_periodic_and_respects_limits() {
    let s = sys(5);
    let (ss, st, _) = s.steady_state(3000).unwrap();
    let (next, _) = s.run_orbit(ss);
    assert!((next.soc - ss.soc).abs() < 1e-9);
    assert!((next.v2 - ss.v2).abs() < 1e-9);
    // CV regulator holds the terminal voltage at or below 4.10 V per cell.
    assert!(st.v_pack_max <= 8.0 * 4.10 + 1e-9);
    assert!(st.i_chg_peak_pack <= 5.0 * 1.6 + 1e-9);
    // Over a periodic orbit charge in equals charge out, so in > out by losses.
    assert!(st.e_in_wh > st.e_out_wh);
    assert!((st.d_chem_wh).abs() < 1e-5);
}

#[test]
fn energy_balance_residual_is_first_order_in_step() {
    let s = sys(5);
    let (ss, _, _) = s.steady_state(3000).unwrap();
    let res = |dt: f64| {
        let mut x = s.clone();
        x.dt = dt;
        x.run_orbit(ss).1.balance_residual()
    };
    let (r2, r1, r05) = (res(2.0), res(1.0), res(0.5));
    assert!(r1 < 1e-4);
    assert!((r2 / r1 - 2.0).abs() < 0.1, "ratio {}", r2 / r1);
    assert!((r1 / r05 - 2.0).abs() < 0.1, "ratio {}", r1 / r05);
}

#[test]
fn sizing_picks_five_strings_and_agrees_with_hand_estimate() {
    let base = sys(1);
    let (np, st) = size_strings(&base, 0.25, 26.0, 20).unwrap();
    assert_eq!(np, 5);
    assert!(st.dod() <= 0.25);
    let (_, st4, _) = sys(4).steady_state(3000).unwrap();
    assert!(st4.dod() > 0.25);
    let q = base.hand_capacity_ah(0.25, 3.70);
    assert_eq!((q / 3.2).ceil() as usize, 5);
    // Hand DoD uses 3.70 V nominal; the pack sits higher in the top quarter
    // of its charge, so the simulated DoD is a little smaller.
    let dod_hand = 0.25 * q / 16.0;
    assert!(st.dod() < dod_hand && st.dod() > 0.9 * dod_hand);
}

#[test]
fn array_power_threshold_brackets_energy_estimate() {
    let s = sys(5);
    let p_min = min_array_power(&s, 0.002).unwrap();
    let mut lo = s.clone();
    lo.p_array = 0.97 * p_min;
    let target = acceptance_limited_soc(&s).unwrap() - 0.002;
    assert!(!recharges(&lo, target));
    // The pure energy balance ignores the voltage-limit taper, so it must
    // come in below the simulated requirement.
    let (_, st, _) = s.steady_state(3000).unwrap();
    let t_e = s.t_eclipse();
    let t_d = s.period() - t_e;
    let p_est = s.p_load + s.p_load * t_e / s.eta_bdr / (st.round_trip() * s.eta_bcr * t_d);
    assert!(p_est < p_min);
}
