// Astron Built this AP
// Project 056 | Analog habitat systems / ECLSS | Multi-zone CO2 ventilation network for a sealed habitat
// Author: Aditya (aditya-astron)
// Date: 2026-11-22 IST
//
// Study driver. `cargo run --release` prints the design study;
// `cargo run --release -- --check` prints key=value lines that
// verify/crosscheck.py recomputes independently with scipy.

use habitat_co2::habitat::*;

const LIMIT_MMHG: f64 = 3.0;
const Q_BASE: f64 = 60.0;
const Q_LEAK: f64 = 3.0;

/// Per-zone (min, max, mean) over one period of the periodic steady state,
/// sampled each minute with the exact map; mean by the trapezoid rule.
fn periodic_stats(h: &Habitat, sched: &[Phase]) -> Vec<(f64, f64, f64)> {
    let n = h.n();
    let mut x = h.periodic_start(sched);
    let mut st: Vec<(f64, f64, f64)> = x.iter().map(|&v| (v, v, 0.0)).collect();
    let mut period = 0.0;
    for ph in sched {
        let steps = (ph.hours * 60.0) as usize;
        let dt = ph.hours / steps as f64;
        let (phi, g) = h.flow_map(&ph.source, dt);
        for _ in 0..steps {
            let nx: Vec<f64> = phi.mul_vec(&x).iter().zip(&g).map(|(u, v)| u + v).collect();
            for k in 0..n {
                st[k].0 = st[k].0.min(nx[k]);
                st[k].1 = st[k].1.max(nx[k]);
                st[k].2 += 0.5 * dt * (x[k] + nx[k]);
            }
            x = nx;
        }
        period += ph.hours;
    }
    st.iter().map(|&(a, b, s)| (a, b, s / period)).collect()
}

/// Peak quarters ppCO2 (mmHg) in the periodic steady state, and when (h after 06:00).
fn quarters_peak(q: f64) -> (f64, f64) {
    let h = reference_habitat(q);
    let sched = reference_schedule();
    let x0 = h.periodic_start(&sched);
    let x_sleep = h.advance(&x0, &sched[0].source, sched[0].hours);
    let (t_day, p_day) = h.phase_peak(&x0, &sched[0], QUARTERS);
    let (t_sl, p_sl) = h.phase_peak(&x_sleep, &sched[1], QUARTERS);
    if p_sl >= p_day { (mmhg(p_sl), sched[0].hours + t_sl) } else { (mmhg(p_day), t_day) }
}

/// Smallest quarters fan flow whose periodic peak stays at the limit.
fn q_min_periodic(limit: f64) -> f64 {
    let (mut lo, mut hi) = (5.0, 500.0);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if quarters_peak(mid).0 > limit { lo = mid } else { hi = mid }
    }
    0.5 * (lo + hi)
}

/// Quarters flow that meets the limit if the sleep sources ran forever
/// (closed-form hub steady state, solved for q).
fn q_min_steady(limit: f64) -> f64 {
    let h = reference_habitat(Q_BASE); // node level does not depend on q
    let c = hub_steady_closed_form(&h, NODE, &sleep_source());
    let s_q = sleep_source()[QUARTERS];
    s_q / (density(limit) - c[NODE])
}

fn hhmm(t_after_0600: f64) -> String {
    let m = ((t_after_0600 * 60.0).round() as i64 + 6 * 60).rem_euclid(24 * 60);
    format!("{:02}:{:02}", m / 60, m % 60)
}

fn main() {
    let check = std::env::args().any(|a| a == "--check");
    let h = reference_habitat(Q_BASE);
    let sched = reference_schedule();
    let n = h.n();
    let total_vol: f64 = h.vol.iter().sum();
    let eff_flow = h.scrub_eff * h.scrub_flow;

    let rates = h.decay_rates();
    let steady_day = h.steady(&day_source());
    let steady_sleep = h.steady(&sleep_source());
    let stats = periodic_stats(&h, &sched);
    let (pk_q, pk_t) = quarters_peak(Q_BASE);
    let q_ss = q_min_steady(LIMIT_MMHG);
    let q_per = q_min_periodic(LIMIT_MMHG);

    // Quarters fan pair fails at lights-out: only a 3 m^3/h door-gap exchange remains.
    let x0 = h.periodic_start(&sched);
    let x_2200 = h.advance(&x0, &sched[0].source, 16.0);
    let failed = reference_habitat(Q_LEAK);
    let t_fan = failed.time_to_limit(&x_2200, &sleep_source(), QUARTERS, density(LIMIT_MMHG), 24.0);

    // Scrubber outage at 06:00 with the day running on.
    let mut no_scrub = reference_habitat(Q_BASE);
    no_scrub.scrub_eff = 0.0;
    let t_scrub: Vec<Option<f64>> = (0..n)
        .map(|k| no_scrub.time_to_limit(&x0, &day_source(), k, density(LIMIT_MMHG), 48.0))
        .collect();

    if check {
        println!("awake_rate={:.12e}", awake_rate());
        for (i, r) in rates.iter().enumerate() {
            println!("rate{}={:.12e}", i, r);
        }
        for k in 0..n {
            println!("steady_day{}={:.12e}", k, steady_day[k]);
            println!("steady_sleep{}={:.12e}", k, steady_sleep[k]);
            println!("periodic0_{}={:.12e}", k, x0[k]);
            println!("periodic2200_{}={:.12e}", k, x_2200[k]);
        }
        println!("t_fan={:.12e}", t_fan.unwrap());
        println!("t_scrub_hab={:.12e}", t_scrub[HAB].unwrap());
        println!("q_min_periodic={:.12e}", q_per);
        return;
    }

    println!("Project 056: CO2 in a five-zone sealed habitat ventilated by fan pairs");
    println!("========================================================================");
    println!("crew CO2: {:.3} kg/d each; awake {:.5} kg/h, asleep {:.5} kg/h; plants take {:.3} kg/h in daylight",
        CO2_PER_DAY, awake_rate(), sleep_rate(), PLANT_UPTAKE);
    println!("ppCO2 per kg/m^3 at {:.2} K: {:.2} mmHg; limit {:.1} mmHg = {:.5} kg/m^3 = {:.0} ppm at 101.325 kPa",
        T_CABIN, mmhg(1.0), LIMIT_MMHG, density(LIMIT_MMHG), LIMIT_MMHG * PA_PER_MMHG / 101325.0 * 1e6);
    println!("scrubber: {:.0} m^3/h at eta {:.2} -> {:.1} m^3/h of CO2-free air in {}",
        h.scrub_flow, h.scrub_eff, eff_flow, h.names[h.scrub_zone]);
    println!();

    println!("1. Zone modes (quarters fan pair {:.0} m^3/h)", Q_BASE);
    let taus: Vec<String> = rates.iter().map(|r| format!("{:.3}", 1.0 / r)).collect();
    println!("   time constants [h]: {}", taus.join("  "));
    println!("   slowest {:.3} h vs. whole-habitat well-mixed V/(eta Qs) = {:.3} h", 1.0 / rates[0], total_vol / eff_flow);
    println!();

    println!("2. Steady states if each phase ran forever [mmHg]");
    println!("   {:<11} {:>8} {:>8}", "zone", "day", "sleep");
    for k in 0..n {
        println!("   {:<11} {:>8.3} {:>8.3}", h.names[k], mmhg(steady_day[k]), mmhg(steady_sleep[k]));
    }
    println!("   well-mixed habitat: day {:.3}, sleep {:.3}",
        mmhg(day_source().iter().sum::<f64>() / eff_flow), mmhg(sleep_source().iter().sum::<f64>() / eff_flow));
    println!();

    println!("3. Periodic steady state over the 24 h schedule [mmHg]");
    println!("   {:<11} {:>7} {:>7} {:>7} {:>8}", "zone", "06:00", "min", "max", "24h mean");
    for k in 0..n {
        println!("   {:<11} {:>7.3} {:>7.3} {:>7.3} {:>8.3}", h.names[k], mmhg(x0[k]), mmhg(stats[k].0), mmhg(stats[k].1), mmhg(stats[k].2));
    }
    let mean_src = (16.0 * day_source().iter().sum::<f64>() + 8.0 * sleep_source().iter().sum::<f64>()) / 24.0;
    println!("   mass balance: eta Qs * mean(c_Hab) = {:.6} kg/h, mean source = {:.6} kg/h",
        eff_flow * stats[HAB].2, mean_src);
    println!("   quarters peak {:.3} mmHg at {}", pk_q, hhmm(pk_t));
    println!();

    println!("4. Sizing the crew-quarters fan pair for a {:.1} mmHg peak", LIMIT_MMHG);
    println!("   {:>8} {:>12} {:>8}", "Q m^3/h", "peak mmHg", "at");
    for q in [20.0, 30.0, 40.0, 50.0, 60.0, 90.0, 120.0] {
        let (p, t) = quarters_peak(q);
        println!("   {:>8.0} {:>12.3} {:>8}", q, p, hhmm(t));
    }
    println!("   sleep steady-state sizing:   Q >= {:.1} m^3/h", q_ss);
    println!("   periodic steady-state sizing: Q >= {:.1} m^3/h  ({:+.0} %)", q_per, 100.0 * (q_per / q_ss - 1.0));
    println!();

    println!("5. Quarters fan pair fails at 22:00 (door gap {:.0} m^3/h left)", Q_LEAK);
    let rise = mmhg(sleep_source()[QUARTERS] / h.vol[QUARTERS]);
    match t_fan {
        Some(t) => println!("   quarters starts at {:.3} mmHg, reaches {:.1} mmHg after {:.0} min (sealed-room rise {:.3} mmHg/h)",
            mmhg(x_2200[QUARTERS]), LIMIT_MMHG, 60.0 * t, rise),
        None => println!("   limit not reached"),
    }
    println!();

    println!("6. Scrubber outage at 06:00, day schedule continuing");
    for k in 0..n {
        match t_scrub[k] {
            Some(t) => println!("   {:<11} reaches {:.1} mmHg after {:.2} h", h.names[k], LIMIT_MMHG, t),
            None => println!("   {:<11} below limit for 48 h", h.names[k]),
        }
    }
    let c_avg: f64 = (0..n).map(|k| h.vol[k] * x0[k]).sum::<f64>() / total_vol;
    println!("   well-mixed estimate from the 06:00 volume-mean level: V (c_lim - c) / S = {:.2} h",
        total_vol * (density(LIMIT_MMHG) - c_avg) / day_source().iter().sum::<f64>());
}
