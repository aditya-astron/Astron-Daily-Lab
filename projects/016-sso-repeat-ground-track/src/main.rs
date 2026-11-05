// Astron Built this AP
// Project 016 | Orbital mechanics | Sun-synchronous repeat ground track design
// Author: Aditya (aditya-astron)
// Date: 2026-11-05 IST

use sso_repeat_track::*;

fn deg(x: f64) -> f64 {
    x.to_degrees()
}

fn main() {
    println!("Sun-synchronous repeat ground track design (J2 secular theory)");
    println!(
        "required node rate = {:.6e} rad/s = {:.5} deg/day\n",
        SUN_RATE,
        deg(SUN_RATE) * 86_400.0
    );

    println!("Published missions (N revs / D days):");
    println!(
        "{:<12} {:>4} {:>3} {:>9} {:>8} {:>9} {:>9} {:>8}  sub-cycle",
        "mission", "N", "D", "a [km]", "h [km]", "i [deg]", "T_N [s]", "dx [km]"
    );
    for (name, n, d) in [("Landsat-8", 233, 16), ("Sentinel-1", 175, 12), ("Sentinel-2", 143, 10)] {
        let o = solve_repeat_sso(n, d).expect("no solution");
        let (sub, east) = o.subcycle();
        println!(
            "{:<12} {:>4} {:>3} {:>9.2} {:>8.2} {:>9.4} {:>9.2} {:>8.1}  {} d, {}",
            name,
            n,
            d,
            o.a,
            o.altitude(),
            deg(o.inc),
            o.nodal_period,
            o.equatorial_spacing(),
            sub,
            if east { "east" } else { "west" }
        );
    }

    println!("\nAll repeat SSOs between 600 and 800 km with cycles up to 7 days:");
    println!(
        "{:>4} {:>3} {:>8} {:>9} {:>9} {:>9}  sub-cycle",
        "N", "D", "h [km]", "i [deg]", "dx [km]", "shift"
    );
    for o in enumerate_repeat_orbits(600.0, 800.0, 7) {
        let (sub, east) = o.subcycle();
        println!(
            "{:>4} {:>3} {:>8.2} {:>9.4} {:>9.1} {:>9.1}  {} d, {}",
            o.revs,
            o.days,
            o.altitude(),
            deg(o.inc),
            o.equatorial_spacing(),
            o.rev_to_rev_shift(),
            sub,
            if east { "east" } else { "west" }
        );
    }

    let o = solve_repeat_sso(233, 16).unwrap();
    println!("\nNumerical J2 propagation of the 233/16 design over one full cycle (RK4, 5 s):");
    let chk = verify_by_propagation(&o, 5.0);
    println!(
        "  nodal period : design {:.4} s, propagated {:.4} s ({:+.2e} rel)",
        o.nodal_period,
        chk.nodal_period,
        chk.nodal_period / o.nodal_period - 1.0
    );
    println!(
        "  node rate    : design {:.5} deg/day, propagated {:.5} deg/day ({:+.2e} rel)",
        deg(SUN_RATE) * 86_400.0,
        deg(chk.raan_rate) * 86_400.0,
        chk.raan_rate / SUN_RATE - 1.0
    );
    println!(
        "  ground track closure after 233 revs: {:+.3} km (track spacing {:.1} km)",
        chk.closure_km,
        o.equatorial_spacing()
    );
}
