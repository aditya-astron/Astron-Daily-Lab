// Astron Built this AP
// Project 035 | Laser/optical communication | Beacon spiral-scan acquisition for a LEO-GEO optical relay link
// Author: Aditya (aditya-astron)
// Date: 2026-11-13 IST

use pat_acquisition::reference::*;
use pat_acquisition::*;

const N_MC: usize = 40_000;

fn urad(x: f64) -> f64 { x * 1e6 }

fn main() {
    let ln = link();
    let pth = ln.p_threshold();
    println!("Beacon acquisition, GEO relay -> LEO user, {:.0} km", ln.range / 1e3);
    println!("  required SNR      {:.3}  (Pfa {:.0e}, Pd {})", ln.required_snr(), ln.pfa, ln.pd);
    println!("  threshold power   {:.3e} W  ({:.2} dBm)", pth, 10.0 * (pth / 1e-3).log10());

    println!("\nBeacon divergence trade (1/e^2 half-angle w):");
    println!("     w [urad]   on-axis margin [dB]   theta_d [urad]");
    let w_opt = ln.optimal_divergence();
    for f in [0.4, 0.6, 0.8, 1.0, 1.25, 1.5, 1.6] {
        let w = f * w_opt;
        let m = 10.0 * (ln.p_rx(w, 0.0) / pth).log10();
        let tag = if f == 1.0 { "  <- optimum" } else { "" };
        println!("   {:9.1}   {:19.2}   {:14.1}{}", urad(w), m, urad(ln.detection_radius(w)), tag);
    }
    let theta_d = ln.detection_radius(w_opt);

    let ideal = design_scan(theta_d, ln.tau, V_MAX);
    let r_scan = scan_radius(SIGMA, P_UNC);
    let dt = ln.tau / 100.0;
    let scan = certify_pitch(&ideal, INNER, ln.tau, r_scan, dt);
    let uniform = certify_pitch(&ideal, 1.0, ln.tau, r_scan, dt);
    println!("\nScan design");
    println!("  pointing sigma    {:.0} urad/axis, cover {} -> R = {:.0} urad", urad(SIGMA), P_UNC, urad(r_scan));
    println!("  dwell-limited     v = {:.3} rad/s  p = {:.0} urad", 2f64.sqrt() * theta_d / ln.tau, urad(2f64.sqrt() * theta_d));
    println!("  straight track    v = {:.3} rad/s  p = {:.0} urad  ({})",
             ideal.speed, urad(ideal.pitch), if ideal.rate_limited { "slew limited" } else { "dwell limited" });
    println!("  certified         p = {:.0} urad ({:.1} % of straight), first turn at {:.0} % speed",
             urad(scan.pitch), 100.0 * scan.pitch / ideal.pitch, 100.0 * INNER);
    println!("  certified, no slow first turn: p = {:.0} urad ({:.1} %)", urad(uniform.pitch), 100.0 * uniform.pitch / ideal.pitch);

    let spiral = Spiral::new(scan.pitch, scan.speed, INNER);
    let track = SpiralTrack::new(&spiral, r_scan, dt);
    let stare = StepStare::new(theta_d, r_scan, ln.tau, T_SETTLE);

    let mc_sp = monte_carlo(SIGMA, N_MC, 35, |x| track.detect(x, theta_d, ln.tau));
    let mc_ss = monte_carlo(SIGMA, N_MC, 35, |x| stare.detect(x));
    let (p_sp, m_sp) = exponential_model(scan.area_rate(), SIGMA, track.duration());
    let m_sp_exact = spiral_mean_time(&spiral, SIGMA, r_scan + scan.pitch);
    let (p_ss, m_ss) = exponential_model(stare.cell_area() / stare.t_spot, SIGMA, stare.duration());

    println!("\nAcquisition, {} Monte Carlo pointing errors", N_MC);
    println!("                     area rate      full scan   P(acq)   E[T] model   E[T] MC   T95 MC");
    println!("  spiral          {:9.3e} rad2/s  {:8.3} s   {:.4}   {:8.3} s   {:7.3} s  {:6.3} s",
             scan.area_rate(), track.duration(), mc_sp.p_acq(), m_sp_exact, mc_sp.mean(), mc_sp.quantile(0.95));
    println!("  step-stare ({:2})  {:9.3e} rad2/s  {:8.3} s   {:.4}   {:8.3} s   {:7.3} s  {:6.3} s",
             stare.spots.len(), stare.cell_area() / stare.t_spot, stare.duration(), mc_ss.p_acq(), m_ss, mc_ss.mean(), mc_ss.quantile(0.95));
    println!("  spiral model: exact arc-length timing with the slow first turn");
    println!("  pure area-rate (exponential) model: spiral E[T] {:.3} s, P {:.4}; step-stare P {:.4}", m_sp, p_sp, p_ss);
    // break-even settle: (3 sqrt3/2) th^2 / (tau + ts) = p v
    let ts_even = stare.cell_area() / scan.area_rate() - ln.tau;
    println!("  step-stare area rate equals the spiral's at a settle time of {:.2} ms", ts_even * 1e3);

    println!("\nBeacon power trade (optimal divergence, same scanner):");
    println!("   P [W]   w [urad]   theta_d [urad]   regime   E[T] [s]   full scan [s]");
    for p in [0.5, 1.0, 2.0, 4.0, 8.0, 16.0] {
        let l = Link { p_beacon: p, ..ln };
        let w = l.optimal_divergence();
        let th = l.detection_radius(w);
        let s = certify_pitch(&design_scan(th, l.tau, V_MAX), INNER, l.tau, r_scan, dt);
        let sp = Spiral::new(s.pitch, s.speed, INNER);
        let t_full = SpiralTrack::new(&sp, r_scan, dt).duration();
        let m = spiral_mean_time(&sp, SIGMA, r_scan + s.pitch);
        println!("   {:5.1}   {:8.0}   {:14.0}   {:6}   {:8.3}   {:13.3}",
                 p, urad(w), urad(th), if s.rate_limited { "slew" } else { "dwell" }, m, t_full);
    }
}
