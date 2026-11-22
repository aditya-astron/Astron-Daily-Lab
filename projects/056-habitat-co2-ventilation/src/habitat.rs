// Astron Built this AP
// Project 056 | Analog habitat systems / ECLSS | Multi-zone CO2 ventilation network for a sealed habitat
// Author: Aditya (aditya-astron)
// Date: 2026-11-22 IST
//
// Well-mixed zone model of CO2 in a habitat ventilated by balanced
// inter-module (IMV) fan pairs and scrubbed in one zone:
//
//   V_i dc_i/dt = sum_j Q_ij (c_j - c_i) - eta Qs c_i [i = scrubber zone] + S_i(t)
//
// c in kg/m^3, V in m^3, Q in m^3/h, S in kg/h, t in h. Written as
// dc/dt = A c + s with A = V^-1 L. Balanced pairs make L symmetric, so A is
// similar to V^-1/2 L V^-1/2 and every mode is real and decaying.
// Sources are piecewise constant over a daily schedule, so each phase is
// propagated exactly with one augmented matrix exponential.

use crate::linalg::{expm, solve, sym_eigenvalues, Mat};

/// Gas constant over CO2 molar mass at cabin temperature, Pa m^3/kg.
pub const T_CABIN: f64 = 295.15;
pub const R_UNIV: f64 = 8.314462618;
pub const M_CO2: f64 = 0.04401;
pub const PA_PER_MMHG: f64 = 133.322368;

/// ppCO2 in mmHg for a CO2 density in kg/m^3 (ideal gas).
pub fn mmhg(c: f64) -> f64 {
    c * R_UNIV * T_CABIN / M_CO2 / PA_PER_MMHG
}

/// CO2 density in kg/m^3 for a partial pressure in mmHg.
pub fn density(pp_mmhg: f64) -> f64 {
    pp_mmhg / mmhg(1.0)
}

#[derive(Clone, Debug)]
pub struct Habitat {
    pub names: Vec<&'static str>,
    pub vol: Vec<f64>,
    /// Balanced fan pairs (i, j, Q): Q m^3/h each way between zones i and j.
    pub links: Vec<(usize, usize, f64)>,
    pub scrub_zone: usize,
    pub scrub_flow: f64,
    pub scrub_eff: f64,
}

#[derive(Clone, Debug)]
pub struct Phase {
    pub name: &'static str,
    pub hours: f64,
    pub source: Vec<f64>, // kg/h per zone, net of plant uptake
}

impl Habitat {
    pub fn n(&self) -> usize {
        self.vol.len()
    }

    /// Exchange-plus-removal matrix L (m^3/h), symmetric.
    pub fn exchange(&self) -> Mat {
        let mut l = Mat::zeros(self.n());
        for &(i, j, q) in &self.links {
            l[(i, j)] += q;
            l[(j, i)] += q;
            l[(i, i)] -= q;
            l[(j, j)] -= q;
        }
        l[(self.scrub_zone, self.scrub_zone)] -= self.scrub_eff * self.scrub_flow;
        l
    }

    /// System matrix A = V^-1 L, 1/h.
    pub fn system(&self) -> Mat {
        let mut a = self.exchange();
        for i in 0..self.n() {
            for j in 0..self.n() {
                a[(i, j)] /= self.vol[i];
            }
        }
        a
    }

    /// Decay rates of the zone modes (1/h, positive, ascending), from the
    /// symmetric form V^-1/2 L V^-1/2.
    pub fn decay_rates(&self) -> Vec<f64> {
        let l = self.exchange();
        let n = self.n();
        let mut b = Mat::zeros(n);
        for i in 0..n {
            for j in 0..n {
                b[(i, j)] = l[(i, j)] / (self.vol[i] * self.vol[j]).sqrt();
            }
        }
        let mut r: Vec<f64> = sym_eigenvalues(&b).into_iter().map(|e| -e).collect();
        r.sort_by(|x, y| x.partial_cmp(y).unwrap());
        r
    }

    pub fn forcing(&self, source: &[f64]) -> Vec<f64> {
        source.iter().zip(&self.vol).map(|(s, v)| s / v).collect()
    }

    /// Steady state for constant sources: A c = -s.
    pub fn steady(&self, source: &[f64]) -> Vec<f64> {
        let rhs: Vec<f64> = self.forcing(source).iter().map(|v| -v).collect();
        solve(&self.system(), &rhs)
    }

    /// Exact affine map over dt hours: c(t+dt) = Phi c(t) + g.
    /// Taken from exp of the augmented generator [[A, s], [0, 0]].
    pub fn flow_map(&self, source: &[f64], dt: f64) -> (Mat, Vec<f64>) {
        let n = self.n();
        let a = self.system();
        let s = self.forcing(source);
        let mut g = Mat::zeros(n + 1);
        for i in 0..n {
            for j in 0..n {
                g[(i, j)] = a[(i, j)] * dt;
            }
            g[(i, n)] = s[i] * dt;
        }
        let e = expm(&g);
        let mut phi = Mat::zeros(n);
        for i in 0..n {
            for j in 0..n {
                phi[(i, j)] = e[(i, j)];
            }
        }
        let gv = (0..n).map(|i| e[(i, n)]).collect();
        (phi, gv)
    }

    pub fn advance(&self, c: &[f64], source: &[f64], dt: f64) -> Vec<f64> {
        let (phi, g) = self.flow_map(source, dt);
        phi.mul_vec(c).iter().zip(&g).map(|(x, y)| x + y).collect()
    }

    /// State at the start of the schedule for the periodic steady state:
    /// one period composes to c_T = M c_0 + b, so c_0 = (I - M)^-1 b.
    /// (I - M) is invertible because every mode decays.
    pub fn periodic_start(&self, schedule: &[Phase]) -> Vec<f64> {
        let n = self.n();
        let mut m = Mat::eye(n);
        let mut b = vec![0.0; n];
        for ph in schedule {
            let (phi, g) = self.flow_map(&ph.source, ph.hours);
            m = phi.mul(&m);
            b = phi.mul_vec(&b).iter().zip(&g).map(|(x, y)| x + y).collect();
        }
        let i_minus_m = Mat::eye(n).add(&m.scale(-1.0));
        solve(&i_minus_m, &b)
    }

    /// Classical RK4 on the ODE, used only as an independent check of the
    /// exponential propagator.
    pub fn rk4(&self, c: &[f64], source: &[f64], dt: f64, steps: usize) -> Vec<f64> {
        let a = self.system();
        let s = self.forcing(source);
        let f = |x: &[f64]| -> Vec<f64> {
            a.mul_vec(x).iter().zip(&s).map(|(u, v)| u + v).collect()
        };
        let axpy = |x: &[f64], k: &[f64], h: f64| -> Vec<f64> {
            x.iter().zip(k).map(|(u, v)| u + h * v).collect()
        };
        let mut x = c.to_vec();
        for _ in 0..steps {
            let k1 = f(&x);
            let k2 = f(&axpy(&x, &k1, dt / 2.0));
            let k3 = f(&axpy(&x, &k2, dt / 2.0));
            let k4 = f(&axpy(&x, &k3, dt));
            for i in 0..x.len() {
                x[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
            }
        }
        x
    }

    /// Maximum of zone k over one phase starting from c0: 1-minute scan to
    /// bracket, then golden-section on the exact propagator.
    pub fn phase_peak(&self, c0: &[f64], ph: &Phase, k: usize) -> (f64, f64) {
        let val = |t: f64| self.advance(c0, &ph.source, t)[k];
        let steps = (ph.hours * 60.0).round() as usize;
        let dt = ph.hours / steps as f64;
        let (phi, g) = self.flow_map(&ph.source, dt);
        let mut x = c0.to_vec();
        let (mut best_i, mut best) = (0usize, x[k]);
        for i in 1..=steps {
            x = phi.mul_vec(&x).iter().zip(&g).map(|(u, v)| u + v).collect();
            if x[k] > best {
                best = x[k];
                best_i = i;
            }
        }
        if best_i == 0 || best_i == steps {
            return (best_i as f64 * dt, best);
        }
        let (mut lo, mut hi) = ((best_i - 1) as f64 * dt, (best_i + 1) as f64 * dt);
        let gr = 0.5 * (5f64.sqrt() - 1.0);
        for _ in 0..60 {
            let m1 = hi - gr * (hi - lo);
            let m2 = lo + gr * (hi - lo);
            if val(m1) > val(m2) {
                hi = m2;
            } else {
                lo = m1;
            }
        }
        let t = 0.5 * (lo + hi);
        (t, val(t))
    }

    /// First time zone k reaches density `limit` from c0 under constant
    /// sources, or None within `horizon` hours. Scan then bisect.
    pub fn time_to_limit(&self, c0: &[f64], source: &[f64], k: usize, limit: f64, horizon: f64) -> Option<f64> {
        if c0[k] >= limit {
            return Some(0.0);
        }
        let dt = 1.0 / 60.0;
        let (phi, g) = self.flow_map(source, dt);
        let mut x = c0.to_vec();
        let mut t = 0.0;
        while t < horizon {
            let nx: Vec<f64> = phi.mul_vec(&x).iter().zip(&g).map(|(u, v)| u + v).collect();
            if nx[k] >= limit {
                let (mut lo, mut hi) = (0.0, dt);
                for _ in 0..60 {
                    let mid = 0.5 * (lo + hi);
                    if self.advance(&x, source, mid)[k] >= limit {
                        hi = mid;
                    } else {
                        lo = mid;
                    }
                }
                return Some(t + 0.5 * (lo + hi));
            }
            x = nx;
            t += dt;
        }
        None
    }
}

// ---- Reference habitat -------------------------------------------------------

pub const HAB: usize = 0;
pub const NODE: usize = 1;
pub const LAB: usize = 2;
pub const GREENHOUSE: usize = 3;
pub const QUARTERS: usize = 4;

/// Metabolic CO2: 1.04 kg/crew/day (NASA BVAD nominal), split so the sleep
/// rate is 0.7 of the awake rate over a 16 h / 8 h day.
pub const CO2_PER_DAY: f64 = 1.04;
pub const SLEEP_FRACTION: f64 = 0.7;
pub fn awake_rate() -> f64 {
    CO2_PER_DAY / (16.0 + 8.0 * SLEEP_FRACTION)
}
pub fn sleep_rate() -> f64 {
    SLEEP_FRACTION * awake_rate()
}

/// Net plant uptake in the greenhouse while lights are on, kg/h.
pub const PLANT_UPTAKE: f64 = 0.02;

/// Five-zone hub layout: every module trades air with the central node.
/// `q_quarters` is the crew-quarters fan pair flow under design.
pub fn reference_habitat(q_quarters: f64) -> Habitat {
    Habitat {
        names: vec!["Hab", "Node", "Lab", "Greenhouse", "Quarters"],
        vol: vec![120.0, 30.0, 80.0, 60.0, 36.0],
        links: vec![(HAB, NODE, 150.0), (LAB, NODE, 100.0), (GREENHOUSE, NODE, 80.0), (QUARTERS, NODE, q_quarters)],
        scrub_zone: HAB,
        scrub_flow: 60.0,
        scrub_eff: 0.85,
    }
}

pub fn day_source() -> Vec<f64> {
    let a = awake_rate();
    vec![2.0 * a, 0.0, a, a - PLANT_UPTAKE, 0.0]
}

pub fn sleep_source() -> Vec<f64> {
    vec![0.0, 0.0, 0.0, 0.0, 4.0 * sleep_rate()]
}

/// 06:00-22:00 awake with greenhouse lights on, 22:00-06:00 all four crew asleep
/// in quarters. t = 0 is 06:00.
pub fn reference_schedule() -> Vec<Phase> {
    vec![
        Phase { name: "day", hours: 16.0, source: day_source() },
        Phase { name: "sleep", hours: 8.0, source: sleep_source() },
    ]
}

/// Steady state of a hub (star) network with the scrubber in a leaf, in
/// closed form. All net CO2 is removed by the scrubber, so
/// c_scrub = S_total / (eta Qs); each link then carries the net source of
/// the subtree behind it, so c rises by (flux / Q) across it.
pub fn hub_steady_closed_form(h: &Habitat, hub: usize, source: &[f64]) -> Vec<f64> {
    let total: f64 = source.iter().sum();
    let mut c = vec![0.0; h.n()];
    let cs = total / (h.scrub_eff * h.scrub_flow);
    c[h.scrub_zone] = cs;
    let q_of = |leaf: usize| -> f64 {
        h.links.iter().find(|&&(i, j, _)| (i == leaf && j == hub) || (j == leaf && i == hub)).unwrap().2
    };
    // Flux from hub to the scrubber leaf is everything not made in that leaf.
    c[hub] = cs + (total - source[h.scrub_zone]) / q_of(h.scrub_zone);
    for k in 0..h.n() {
        if k != hub && k != h.scrub_zone {
            c[k] = c[hub] + source[k] / q_of(k);
        }
    }
    c
}
