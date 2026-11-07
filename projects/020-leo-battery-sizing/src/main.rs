// Astron Built this AP
// Project 020 | Spacecraft power systems | LEO Li-ion battery sizing with a 2RC equivalent-circuit model
// Author: Aditya (aditya-astron)
// Date: 2026-11-07 IST

use leo_battery::*;

fn baseline() -> PowerSystem {
    PowerSystem {
        cell: Cell::nmc_18650(),
        ns: 8,
        np: 1,
        p_load: 180.0,
        p_array: 400.0,
        eta_bcr: 0.95,
        eta_bdr: 0.95,
        v_cell_max: 4.10,
        i_chg_max: 0.5 * 3.2, // C/2 per cell
        alt_m: 500e3,
        beta_rad: 0.0,
        dt: 1.0,
    }
}

fn main() {
    let dod_max = 0.25;
    let v_bus_floor = 26.0; // BDR input floor, 8s pack
    let v_cell_nom = 3.70;

    let base = baseline();
    let t_orb = base.period();
    println!("=== Orbit: 500 km circular ===");
    println!("period                 {:8.1} s  ({:.2} min)", t_orb, t_orb / 60.0);
    println!("no-eclipse beta        {:8.2} deg", beta_no_eclipse(base.alt_m).to_degrees());
    println!("  beta[deg]  eclipse[min]  fraction");
    for b in [0.0, 20.0, 40.0, 60.0, 65.0, 70.0] {
        let f = eclipse_fraction(base.alt_m, f64::to_radians(b));
        println!("  {:8.1}  {:12.2}  {:8.4}", b, f * t_orb / 60.0, f);
    }

    println!("\n=== Sizing (beta = 0, load {:.0} W, DoD <= {:.0}%) ===", base.p_load, dod_max * 100.0);
    let q_hand = base.hand_capacity_ah(dod_max, v_cell_nom);
    let np_hand = (q_hand / base.cell.q_ah).ceil() as usize;
    println!("eclipse battery energy {:8.2} Wh", base.p_load * base.t_eclipse() / 3600.0 / base.eta_bdr);
    println!("hand capacity          {:8.2} Ah  -> {} strings of {:.1} Ah", q_hand, np_hand, base.cell.q_ah);

    let e_ecl_wh = base.p_load * base.t_eclipse() / 3600.0 / base.eta_bdr;
    println!("  np  hand DoD[%]  DoD[%]  Vmin[V]  Vmax[V]  Ichg[A]  t_CV[min]  round-trip");
    for np in 1..=7 {
        let mut s = base.clone();
        s.np = np;
        match s.steady_state(3000) {
            Some((_, st, _)) => println!(
                "  {:2}  {:11.2}  {:6.2}  {:7.2}  {:7.2}  {:7.2}  {:9.1}  {:9.4}",
                np,
                100.0 * e_ecl_wh / (base.ns as f64 * v_cell_nom * base.cell.q_ah * np as f64),
                st.dod() * 100.0,
                st.v_pack_min,
                st.v_pack_max,
                st.i_chg_peak_pack,
                st.t_cv_s / 60.0,
                st.round_trip()
            ),
            None => println!("  {:2}  battery runs down: eclipse draw exceeds capacity or C/2 recharge", np),
        }
    }
    let (np, st) = size_strings(&base, dod_max, v_bus_floor, 20).expect("no feasible pack");
    println!("simulated sizing       {} strings ({}s{}p, {:.1} Ah, {:.0} Wh nameplate)",
        np, base.ns, np, np as f64 * base.cell.q_ah,
        np as f64 * base.cell.q_ah * base.ns as f64 * v_cell_nom);

    let mut sys = base.clone();
    sys.np = np;
    let (s_ss, _, n_orb) = sys.steady_state(3000).unwrap();
    println!("\n=== Steady-state orbit, {}s{}p, array {:.0} W ===", sys.ns, np, sys.p_array);
    println!("converged after        {:8} orbits", n_orb);
    println!("eclipse-entry SOC      {:8.4}  (CV-full {:.4})", s_ss.soc, sys.soc_full());
    println!("DoD                    {:8.2} %", st.dod() * 100.0);
    println!("pack voltage           {:8.2} .. {:.2} V", st.v_pack_min, st.v_pack_max);
    println!("peak discharge/charge  {:8.2} / {:.2} A", st.i_dis_peak_pack, st.i_chg_peak_pack);
    println!("time on voltage limit  {:8.1} min", st.t_cv_s / 60.0);
    println!("energy out / in        {:8.2} / {:.2} Wh", st.e_out_wh, st.e_in_wh);
    println!("Joule loss             {:8.2} Wh", st.e_loss_wh);
    println!("round-trip efficiency  {:8.4}", st.round_trip());
    println!("array energy shunted   {:8.2} Wh", st.e_shunt_wh);

    println!("\n  energy-balance residual vs step size (first order expected)");
    for dt in [4.0, 2.0, 1.0, 0.5] {
        let mut s = sys.clone();
        s.dt = dt;
        let (_, so) = s.run_orbit(s_ss);
        println!("    dt = {:4.1} s   residual {:.3e}", dt, so.balance_residual());
    }

    println!("\n=== Array power needed to recharge every orbit ===");
    let soc_acc = acceptance_limited_soc(&sys).unwrap();
    let p_min = min_array_power(&sys, 0.002).unwrap();
    let t_e = sys.t_eclipse();
    let t_d = sys.period() - t_e;
    // Energy-balance estimate: the sunlit surplus must refill what the
    // eclipse took, through the BCR and the battery's own round-trip loss.
    let e_batt_out = sys.p_load * t_e / sys.eta_bdr;
    let p_est = sys.p_load + e_batt_out / (st.round_trip() * sys.eta_bcr * t_d);
    println!("acceptance-limited entry SOC        {:7.4}", soc_acc);
    println!("bisected (entry SOC within 0.2%)    {:7.1} W", p_min);
    println!("energy-balance estimate             {:7.1} W", p_est);
    println!("ratio                               {:7.3}", p_min / p_est);

    println!("\n=== {}s{}p across beta (array {:.0} W) ===", sys.ns, np, sys.p_array);
    println!("  beta[deg]  DoD[%]  Vmin[V]  t_CV[min]");
    for b in [0.0, 30.0, 50.0, 60.0, 65.0] {
        let mut s = sys.clone();
        s.beta_rad = f64::to_radians(b);
        let (_, so, _) = s.steady_state(3000).unwrap();
        println!("  {:8.1}  {:6.2}  {:7.2}  {:9.1}", b, so.dod() * 100.0, so.v_pack_min, so.t_cv_s / 60.0);
    }
}
