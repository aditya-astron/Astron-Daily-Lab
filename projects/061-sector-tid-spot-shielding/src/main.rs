// Astron Built this AP
// Project 061 | Radiation shielding | Sector-based TID analysis and spot-shield sizing for an avionics box
// Author: Aditya (aditya-astron)
// Date: 2026-11-23 IST
//
// Design study: ray-traced sector TID for four parts in a GEO avionics box,
// checks of the sector integrator against closed forms, and the mass cost
// of four ways to bring the failing parts inside a 2x radiation design
// margin (RDM): spot-shield cans, cans with a backside plate, a +X wall
// doubler, and uniform thickening of the box walls.

use sector_tid::dose::*;
use sector_tid::geometry::Scene;
use sector_tid::model::*;
use sector_tid::sector::*;

const RDM: f64 = 2.0;
const N_DIRS: usize = 200_000;

fn scene_with(extra: f64, doubler: f64, cans: &[Scene]) -> Scene {
    let mut s = bus_scene();
    s.extend(&avionics_box(extra, doubler));
    for c in cans {
        s.extend(c);
    }
    s
}

fn main() {
    let dd = DoseDepth::geo_15yr_illustrative();
    let dirs = fibonacci_dirs(N_DIRS);
    let parts = parts();

    println!("1. Dose-depth curve (solid Al sphere, illustrative 15-yr GEO)");
    println!("   mm Al   g/cm^2   krad(Si)");
    for mm in [0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0] {
        let t = mm_al_to_areal(mm);
        println!("   {:5.1}   {:6.3}   {:9.2}", mm, t, dd.eval(t));
    }

    println!("\n2. Sector integrator against closed forms (N = {N_DIRS} Fibonacci directions)");
    {
        let mut s = Scene::new();
        s.add_sphere_shell([0.0; 3], 20.0, 0.3, AL_RHO);
        let got = sector_dose(&s, &[0.0; 3], &dd, &dirs);
        let want = dd.eval(0.3 * AL_RHO);
        println!("   sphere shell, centre       {:10.6}  exact {:10.6}  rel {:+.1e}", got, want, got / want - 1.0);

        let got = sector_dose(&s, &[0.0, 0.0, 12.0], &dd, &dirs);
        let want = offset_shell_reference(&dd, 19.7, 20.0, AL_RHO, 12.0, 1e-12);
        println!("   sphere shell, d = 0.6 R    {:10.6}  quad  {:10.6}  rel {:+.1e}", got, want, got / want - 1.0);

        let w = 0.54;
        let mut s = Scene::new();
        s.add_box([-1e6, -1e6, 1.0], [1e6, 1e6, 1.0 + w / AL_RHO], AL_RHO);
        s.add_box([-1e6, -1e6, -1.0 - w / AL_RHO], [1e6, 1e6, -1.0], AL_RHO);
        let got = sector_dose(&s, &[0.0; 3], &dd, &dirs);
        let want = dd.two_sided_slab(w);
        println!("   two infinite plates        {:10.6}  E2    {:10.6}  rel {:+.1e}", got, want, got / want - 1.0);

        let mut s = Scene::new();
        s.add_box_shell([-10.0; 3], [10.0; 3], 0.2, AL_RHO);
        let got = sector_dose(&s, &[0.0; 3], &dd, &dirs);
        let want = cube_centre_reference(&dd, 0.2 * AL_RHO, 1e-13);
        println!("   cubic shell, centre        {:10.6}  quad  {:10.6}  rel {:+.1e}", got, want, got / want - 1.0);
    }

    let base = scene_with(0.0, 0.0, &[]);
    let p2 = parts[1].point();
    println!("\n3. Convergence at {} in the full model", parts[1].name);
    let t_ref = thickness_map(&base, &p2, &fibonacci_dirs(4_000_000));
    let d_ref = fold(&dd, &t_ref);
    println!("   reference (4M Fibonacci) {:.6} krad", d_ref);
    println!("        N    Fibonacci err    random err   random 1-sigma");
    for n in [1_000, 10_000, 100_000, 1_000_000] {
        let fib = sector_dose(&base, &p2, &dd, &fibonacci_dirs(n));
        let t = thickness_map(&base, &p2, &random_dirs(n, 12345));
        let doses: Vec<f64> = t.iter().map(|&x| dd.eval(x)).collect();
        let mean = doses.iter().sum::<f64>() / n as f64;
        let var = doses.iter().map(|d| (d - mean) * (d - mean)).sum::<f64>() / (n - 1) as f64;
        println!(
            "   {:7}   {:+12.2e}   {:+11.2e}   {:12.2e}",
            n,
            fib / d_ref - 1.0,
            mean / d_ref - 1.0,
            (var / n as f64).sqrt() / d_ref
        );
    }

    println!("\n4. Baseline sector doses, requirement D * {RDM} <= capability");
    println!("   part                       D krad  eq.sphere  t_min  D(t_min)  t<1 g/cm^2   cap  margin");
    let mut doses0 = Vec::new();
    for p in &parts {
        let t = thickness_map(&base, &p.point(), &dirs);
        let d = fold(&dd, &t);
        doses0.push(d);
        let tmin = percentile(&t, 0.0);
        let thin = t.iter().filter(|&&x| x < 1.0).count() as f64 / t.len() as f64;
        println!(
            "   {}  {:6.2}  {:5.2} mm  {:5.3}  {:8.1}      {:5.1} %  {:4.0}  {:5.2} {}",
            p.name,
            d,
            areal_to_mm_al(dd.equivalent_thickness(d)),
            tmin,
            dd.eval(tmin),
            100.0 * thin,
            p.capability,
            p.capability / d,
            if p.capability / d >= RDM { "pass" } else { "FAIL" }
        );
    }

    println!("\n5. Where the dose comes from (share of D by avionics-box face crossed)");
    let names = ["+X (panel)", "-X (tank)", "+Y", "-Y", "+Z", "-Z"];
    for (i, p) in parts.iter().enumerate() {
        if i == 0 || i == 2 {
            continue;
        }
        let pt = p.point();
        let mut share = [0.0; 6];
        let mut omega = [0.0; 6];
        for u in &dirs {
            let f = exit_face(&pt, u);
            share[f] += dd.eval(base.areal_density(&pt, u));
            omega[f] += 1.0;
        }
        let tot: f64 = share.iter().sum();
        print!("   {}:", &p.name[..7]);
        for f in 0..6 {
            print!("  {} {:4.1}% ({:4.1}% of 4pi)", names[f], 100.0 * share[f] / tot, 100.0 * omega[f] / N_DIRS as f64);
            if f == 2 {
                print!("\n           ");
            }
        }
        println!();
    }

    println!("\n6. Mitigation trade for the failing parts (Al added, 1 um bisection)");
    let failing: Vec<usize> = (0..parts.len()).filter(|&i| parts[i].capability / doses0[i] < RDM).collect();
    let target = |i: usize| parts[i].capability / RDM;

    // a, b: spot cans sized per part (each sized with only its own can).
    let mut can_designs: Vec<(bool, Vec<(usize, f64)>, f64)> = Vec::new();
    for backside in [false, true] {
        let mut ws = Vec::new();
        let mut mass = 0.0;
        for &i in &failing {
            let p = &parts[i];
            let pt = p.point();
            let f = |w: f64| sector_dose(&scene_with(0.0, 0.0, &[spot_can(p, w, backside)]), &pt, &dd, &dirs);
            match bisect_decreasing(f, target(i), 0.0, 1.0, 1e-4) {
                Some(w) => {
                    let m = spot_can(p, w, backside).mass();
                    println!(
                        "   {} can{} on {}: wall {:.2} mm, {:6.1} g",
                        if backside { "6-sided" } else { "5-sided" },
                        if backside { " + back plate" } else { "" },
                        &p.name[..7],
                        10.0 * w,
                        m
                    );
                    ws.push((i, w));
                    mass += m;
                }
                None => println!("   {} can on {}: not reachable with 10 mm walls", if backside { "6-sided" } else { "5-sided" }, &p.name[..7]),
            }
        }
        can_designs.push((backside, ws, mass));
    }

    // c: +X doubler.
    let box0 = avionics_box(0.0, 0.0).mass();
    let dbl = |d: f64, i: usize| sector_dose(&scene_with(0.0, d, &[]), &parts[i].point(), &dd, &dirs);
    let mut d_need: f64 = 0.0;
    let mut dbl_ok = true;
    for &i in &failing {
        match bisect_decreasing(|d| dbl(d, i), target(i), 0.0, 3.0, 1e-4) {
            Some(d) => {
                println!("   +X doubler for {}: {:.2} mm", &parts[i].name[..7], 10.0 * d);
                d_need = d_need.max(d);
            }
            None => {
                println!(
                    "   +X doubler for {}: not enough even at 30 mm ({:.2} krad, needs {:.2})",
                    &parts[i].name[..7],
                    dbl(3.0, i),
                    target(i)
                );
                dbl_ok = false;
            }
        }
    }
    println!(
        "   +X doubler {:.2} mm: {:.1} g{}",
        10.0 * d_need,
        avionics_box(0.0, d_need).mass() - box0,
        if dbl_ok { "" } else { " (does not close the design)" }
    );

    // d: uniform wall thickening.
    let thick = |e: f64, i: usize| sector_dose(&scene_with(e, 0.0, &[]), &parts[i].point(), &dd, &dirs);
    let mut e_need: f64 = 0.0;
    for &i in &failing {
        let e = bisect_decreasing(|e| thick(e, i), target(i), 0.0, 0.8, 1e-4).expect("within 8 mm of card clearance");
        println!("   uniform wall for {}: +{:.2} mm", &parts[i].name[..7], 10.0 * e);
        e_need = e_need.max(e);
    }
    let m_thick = avionics_box(e_need, 0.0).mass() - box0;
    println!("   uniform wall {:.2} -> {:.2} mm: {:.1} g", 10.0 * BOX_WALL, 10.0 * (BOX_WALL + e_need), m_thick);

    println!("\n7. Chosen design: 5-sided cans, all parts re-evaluated together");
    let (_, ws, mass) = &can_designs[0];
    let cans: Vec<Scene> = ws.iter().map(|&(i, w)| spot_can(&parts[i], w, false)).collect();
    let fin = scene_with(0.0, 0.0, &cans);
    println!("   part                       before   after  margin");
    for (i, p) in parts.iter().enumerate() {
        let d = sector_dose(&fin, &p.point(), &dd, &dirs);
        println!("   {}  {:6.2}  {:6.2}  {:5.2} {}", p.name, doses0[i], d, p.capability / d, if p.capability / d >= RDM - 1e-6 { "pass" } else { "FAIL" });
    }
    println!("   added mass {:.1} g vs {:.1} g for uniform walls ({:.0}x)", mass, m_thick, m_thick / mass);
}
