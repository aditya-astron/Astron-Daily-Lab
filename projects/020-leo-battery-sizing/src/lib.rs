// Astron Built this AP
// Project 020 | Spacecraft power systems | LEO Li-ion battery sizing with a 2RC equivalent-circuit model
// Author: Aditya (aditya-astron)
// Date: 2026-11-07 IST
//
// Battery model: Thevenin equivalent circuit with two RC branches per cell,
//
//     V_t = OCV(SOC) - V1 - V2 - I R0          (I > 0 is discharge)
//     dVk/dt = -Vk / (Rk Ck) + I / Ck
//     dSOC/dt = -I / Q
//
// The pack is ns cells in series by np strings in parallel, all cells assumed
// matched, so the simulation runs on one cell and scales power by ns*np.
// Current is held constant over a step (the regulator sets it at the start),
// so the RC branches are discretised exactly with a zero-order hold.

pub const MU_EARTH: f64 = 3.986_004_418e14; // m^3/s^2
pub const R_EARTH: f64 = 6_378_137.0; // m

// ---------------------------------------------------------------- orbit

/// Circular orbit period [s].
pub fn orbit_period(alt_m: f64) -> f64 {
    let a = R_EARTH + alt_m;
    2.0 * std::f64::consts::PI * (a * a * a / MU_EARTH).sqrt()
}

/// Fraction of a circular orbit spent in Earth's cylindrical shadow.
///
/// The orbit dips into the shadow cylinder when the projection of the
/// position on the plane normal to the Sun is shorter than R_E. For a Sun
/// elevation beta above the orbit plane that gives
///     cos(pi f) = sqrt(h^2 + 2 R h) / (r cos beta),
/// and no eclipse once sin(beta) >= R / r.
pub fn eclipse_fraction(alt_m: f64, beta_rad: f64) -> f64 {
    let r = R_EARTH + alt_m;
    let x = (alt_m * alt_m + 2.0 * R_EARTH * alt_m).sqrt() / (r * beta_rad.cos());
    if x >= 1.0 {
        0.0
    } else {
        x.acos() / std::f64::consts::PI
    }
}

/// Beta angle above which a circular orbit has no eclipse [rad].
pub fn beta_no_eclipse(alt_m: f64) -> f64 {
    (R_EARTH / (R_EARTH + alt_m)).asin()
}

// ---------------------------------------------------------------- cell

/// Generic 18650-class NMC cell. OCV table is a typical low-rate curve.
#[derive(Clone, Debug)]
pub struct Cell {
    pub q_ah: f64,
    pub r0: f64,
    pub r1: f64,
    pub c1: f64,
    pub r2: f64,
    pub c2: f64,
    pub soc_pts: Vec<f64>,
    pub ocv_pts: Vec<f64>,
}

impl Cell {
    pub fn nmc_18650() -> Self {
        Cell {
            q_ah: 3.2,
            r0: 0.030,
            r1: 0.015,
            c1: 20.0 / 0.015, // tau1 = 20 s, charge transfer / double layer
            r2: 0.020,
            c2: 600.0 / 0.020, // tau2 = 600 s, solid-phase diffusion
            soc_pts: vec![0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
            ocv_pts: vec![3.00, 3.35, 3.46, 3.56, 3.62, 3.67, 3.73, 3.82, 3.91, 3.99, 4.08, 4.18],
        }
    }

    pub fn q_as(&self) -> f64 {
        self.q_ah * 3600.0
    }

    fn seg(&self, soc: f64) -> usize {
        let n = self.soc_pts.len();
        let mut k = 0;
        while k < n - 2 && soc > self.soc_pts[k + 1] {
            k += 1;
        }
        k
    }

    /// Open-circuit voltage, linear interpolation (linear extrapolation outside).
    pub fn ocv(&self, soc: f64) -> f64 {
        let k = self.seg(soc);
        let (s0, s1) = (self.soc_pts[k], self.soc_pts[k + 1]);
        let (v0, v1) = (self.ocv_pts[k], self.ocv_pts[k + 1]);
        v0 + (v1 - v0) * (soc - s0) / (s1 - s0)
    }

    /// W(soc) = integral_0^soc OCV(s) ds, exact for the piecewise-linear table.
    /// Chemical energy stored per cell is Q * W(soc).
    pub fn ocv_integral(&self, soc: f64) -> f64 {
        let mut w = 0.0;
        let k_end = self.seg(soc);
        for k in 0..k_end {
            let ds = self.soc_pts[k + 1] - self.soc_pts[k];
            w += 0.5 * (self.ocv_pts[k] + self.ocv_pts[k + 1]) * ds;
        }
        let s0 = self.soc_pts[k_end];
        w + 0.5 * (self.ocv_pts[k_end] + self.ocv(soc)) * (soc - s0)
    }

    /// SOC at which OCV equals v (bisection; OCV is monotone).
    pub fn soc_at_ocv(&self, v: f64) -> f64 {
        let (mut lo, mut hi) = (0.0, 1.0);
        for _ in 0..80 {
            let m = 0.5 * (lo + hi);
            if self.ocv(m) < v {
                lo = m
            } else {
                hi = m
            }
        }
        0.5 * (lo + hi)
    }
}

/// Cell electrical state.
#[derive(Clone, Copy, Debug)]
pub struct CellState {
    pub soc: f64,
    pub v1: f64,
    pub v2: f64,
}

impl CellState {
    pub fn rested(soc: f64) -> Self {
        CellState { soc, v1: 0.0, v2: 0.0 }
    }
}

pub fn terminal_voltage(c: &Cell, s: &CellState, i: f64) -> f64 {
    c.ocv(s.soc) - s.v1 - s.v2 - i * c.r0
}

/// Discharge current that delivers power p at the terminals:
///     p = (E - I R0) I,  E = OCV - V1 - V2
/// Smaller root of the quadratic (the stable operating point). None if
/// p exceeds the maximum transferable power E^2 / (4 R0).
pub fn current_for_discharge_power(c: &Cell, s: &CellState, p: f64) -> Option<f64> {
    let e = c.ocv(s.soc) - s.v1 - s.v2;
    let disc = e * e - 4.0 * c.r0 * p;
    if disc < 0.0 {
        return None;
    }
    // Rationalised form avoids cancellation for small p.
    Some(2.0 * p / (e + disc.sqrt()))
}

/// Charge current magnitude that absorbs power p: p = (E + Ic R0) Ic.
pub fn current_for_charge_power(c: &Cell, s: &CellState, p: f64) -> f64 {
    let e = c.ocv(s.soc) - s.v1 - s.v2;
    2.0 * p / (e + (e * e + 4.0 * c.r0 * p).sqrt())
}

/// Advance the cell by dt at constant current i (exact ZOH for the RC states).
pub fn step_cell(c: &Cell, s: &CellState, i: f64, dt: f64) -> CellState {
    let a1 = (-dt / (c.r1 * c.c1)).exp();
    let a2 = (-dt / (c.r2 * c.c2)).exp();
    CellState {
        soc: s.soc - i * dt / c.q_as(),
        v1: a1 * s.v1 + c.r1 * (1.0 - a1) * i,
        v2: a2 * s.v2 + c.r2 * (1.0 - a2) * i,
    }
}

/// Exact integral of Vk^2 / Rk over a ZOH step: Joule loss in an RC branch.
fn rc_loss(v0: f64, i: f64, r: f64, cap: f64, dt: f64) -> f64 {
    let tau = r * cap;
    let vs = i * r;
    let d = v0 - vs;
    let a = (-dt / tau).exp();
    (vs * vs * dt + 2.0 * vs * d * tau * (1.0 - a) + d * d * 0.5 * tau * (1.0 - a * a)) / r
}

// ---------------------------------------------------------------- system

#[derive(Clone, Debug)]
pub struct PowerSystem {
    pub cell: Cell,
    pub ns: usize,
    pub np: usize,
    pub p_load: f64,     // regulated bus load [W]
    pub p_array: f64,    // array power available in sunlight [W]
    pub eta_bcr: f64,    // battery charge regulator efficiency
    pub eta_bdr: f64,    // battery discharge regulator efficiency
    pub v_cell_max: f64, // charge voltage limit per cell (CV setpoint)
    pub i_chg_max: f64,  // charge current limit per cell [A]
    pub alt_m: f64,
    pub beta_rad: f64,
    pub dt: f64,
}

#[derive(Clone, Debug, Default)]
pub struct OrbitStats {
    pub soc_start: f64,
    pub soc_min: f64,
    pub soc_end: f64,
    pub v_pack_min: f64,
    pub v_pack_max: f64,
    pub i_chg_peak_pack: f64,
    pub i_dis_peak_pack: f64,
    pub t_cv_s: f64,          // time the voltage limit sets the charge current
    pub e_out_wh: f64,        // battery terminal energy delivered
    pub e_in_wh: f64,         // battery terminal energy accepted
    pub e_loss_wh: f64,       // R0 + RC Joule losses
    pub d_chem_wh: f64,       // change in stored chemical energy
    pub d_cap_wh: f64,        // change in RC capacitor energy
    pub e_shunt_wh: f64,      // array power left unused
    pub collapsed: bool,      // load exceeded the deliverable power
}

impl OrbitStats {
    pub fn dod(&self) -> f64 {
        self.soc_start - self.soc_min
    }
    /// Battery round-trip energy efficiency over the orbit (terminal out / in).
    pub fn round_trip(&self) -> f64 {
        self.e_out_wh / self.e_in_wh
    }
    /// Energy balance residual, normalised by energy delivered.
    /// in - out = d_chem + d_cap + loss.
    pub fn balance_residual(&self) -> f64 {
        let lhs = self.e_in_wh - self.e_out_wh;
        let rhs = self.d_chem_wh + self.d_cap_wh + self.e_loss_wh;
        (lhs - rhs).abs() / self.e_out_wh
    }
}

impl PowerSystem {
    pub fn period(&self) -> f64 {
        orbit_period(self.alt_m)
    }
    pub fn t_eclipse(&self) -> f64 {
        eclipse_fraction(self.alt_m, self.beta_rad) * self.period()
    }
    pub fn n_cells(&self) -> f64 {
        (self.ns * self.np) as f64
    }
    pub fn soc_full(&self) -> f64 {
        self.cell.soc_at_ocv(self.v_cell_max)
    }

    fn cap_energy(&self, s: &CellState) -> f64 {
        0.5 * self.cell.c1 * s.v1 * s.v1 + 0.5 * self.cell.c2 * s.v2 * s.v2
    }

    /// Simulate one orbit starting at eclipse entry.
    pub fn run_orbit(&self, s0: CellState) -> (CellState, OrbitStats) {
        let c = &self.cell;
        let t_orb = self.period();
        let t_ecl = self.t_eclipse();
        let n = (t_orb / self.dt).round() as usize;
        let dt = t_orb / n as f64;
        let k = self.n_cells();
        let to_wh = k / 3600.0;

        let mut st = OrbitStats {
            soc_start: s0.soc,
            soc_min: s0.soc,
            v_pack_min: f64::INFINITY,
            v_pack_max: f64::NEG_INFINITY,
            ..Default::default()
        };
        let mut s = s0;
        let e_chem0 = c.q_as() * c.ocv_integral(s.soc);
        let e_cap0 = self.cap_energy(&s);
        let (mut e_out, mut e_in, mut e_loss, mut e_shunt) = (0.0, 0.0, 0.0, 0.0);

        for j in 0..n {
            let t = (j as f64 + 0.5) * dt;
            let sunlit = t >= t_ecl;
            let p_sa = if sunlit { self.p_array } else { 0.0 };
            let surplus = p_sa - self.p_load;

            let mut cv = false;
            let i = if surplus < 0.0 {
                let p_cell = -surplus / self.eta_bdr / k;
                match current_for_discharge_power(c, &s, p_cell) {
                    Some(i) => i,
                    None => {
                        st.collapsed = true;
                        break;
                    }
                }
            } else {
                let p_cell = surplus * self.eta_bcr / k;
                let i_pow = current_for_charge_power(c, &s, p_cell);
                let e = c.ocv(s.soc) - s.v1 - s.v2;
                let i_volt = ((self.v_cell_max - e) / c.r0).max(0.0);
                let mut ic = i_pow.min(self.i_chg_max);
                if i_volt < ic {
                    ic = i_volt;
                    cv = true;
                }
                // Array power not taken by the battery goes to the shunt.
                let p_taken = (e + ic * c.r0) * ic * k / self.eta_bcr;
                e_shunt += (surplus - p_taken).max(0.0) * dt;
                -ic
            };

            let vt = terminal_voltage(c, &s, i);
            let vp = vt * self.ns as f64;
            st.v_pack_min = st.v_pack_min.min(vp);
            st.v_pack_max = st.v_pack_max.max(vp);
            let ip = i * self.np as f64;
            st.i_dis_peak_pack = st.i_dis_peak_pack.max(ip);
            st.i_chg_peak_pack = st.i_chg_peak_pack.max(-ip);
            if cv {
                st.t_cv_s += dt;
            }

            // Terminal energy as a power monitor would integrate it:
            // sampled at the start of the step, held for dt (first order in dt).
            let p_term = vt * i;
            if p_term > 0.0 {
                e_out += p_term * dt;
            } else {
                e_in += -p_term * dt;
            }
            e_loss += i * i * c.r0 * dt
                + rc_loss(s.v1, i, c.r1, c.c1, dt)
                + rc_loss(s.v2, i, c.r2, c.c2, dt);

            s = step_cell(c, &s, i, dt);
            st.soc_min = st.soc_min.min(s.soc);
        }

        st.soc_end = s.soc;
        st.e_out_wh = e_out * to_wh;
        st.e_in_wh = e_in * to_wh;
        st.e_loss_wh = e_loss * to_wh;
        st.d_chem_wh = (c.q_as() * c.ocv_integral(s.soc) - e_chem0) * to_wh;
        st.d_cap_wh = (self.cap_energy(&s) - e_cap0) * to_wh;
        st.e_shunt_wh = e_shunt / 3600.0;
        (s, st)
    }

    /// Repeat orbits from a full, rested battery until the eclipse-entry state
    /// repeats (periodic steady state). Returns None on collapse or if SOC is
    /// still drifting after max_orbits (energy-negative orbit).
    pub fn steady_state(&self, max_orbits: usize) -> Option<(CellState, OrbitStats, usize)> {
        let mut s = CellState::rested(self.soc_full());
        for n in 1..=max_orbits {
            let (s_next, st) = self.run_orbit(s);
            if st.collapsed || s_next.soc <= 0.0 {
                return None;
            }
            let ds = (s_next.soc - s.soc).abs() + (s_next.v1 - s.v1).abs() + (s_next.v2 - s.v2).abs();
            s = s_next;
            if ds < 1e-10 {
                // One more orbit from the converged state for clean statistics.
                let (_, st) = self.run_orbit(s);
                return Some((s, st, n));
            }
        }
        None
    }

    /// Hand sizing: Ah needed so that the eclipse draw is dod_max of capacity,
    /// at nominal pack voltage.
    pub fn hand_capacity_ah(&self, dod_max: f64, v_cell_nom: f64) -> f64 {
        let e_wh = self.p_load * self.t_eclipse() / 3600.0 / self.eta_bdr;
        e_wh / (self.ns as f64 * v_cell_nom * dod_max)
    }
}

/// Smallest number of parallel strings meeting the DoD and minimum-voltage
/// limits in periodic steady state.
pub fn size_strings(base: &PowerSystem, dod_max: f64, v_pack_min: f64, np_max: usize) -> Option<(usize, OrbitStats)> {
    for np in 1..=np_max {
        let mut sys = base.clone();
        sys.np = np;
        if let Some((_, st, _)) = sys.steady_state(3000) {
            if st.dod() <= dod_max && st.v_pack_min >= v_pack_min {
                return Some((np, st));
            }
        }
    }
    None
}

/// Does the battery reach at least soc_target by eclipse entry, in periodic
/// steady state?
pub fn recharges(sys: &PowerSystem, soc_target: f64) -> bool {
    match sys.steady_state(3000) {
        Some((s, _, _)) => s.soc >= soc_target,
        None => false,
    }
}

/// Eclipse-entry SOC when array power is effectively unlimited, so charging
/// is set only by the current limit and the voltage-limit taper.
pub fn acceptance_limited_soc(base: &PowerSystem) -> Option<f64> {
    let mut sys = base.clone();
    sys.p_array = 20.0 * base.p_load;
    sys.steady_state(3000).map(|(s, _, _)| s.soc)
}

/// Minimum array power whose steady-state eclipse-entry SOC is within tol of
/// the acceptance-limited value: beyond it, more array buys nothing (bisection).
pub fn min_array_power(base: &PowerSystem, tol: f64) -> Option<f64> {
    let target = acceptance_limited_soc(base)? - tol;
    let mut sys = base.clone();
    let (mut lo, mut hi) = (base.p_load, 20.0 * base.p_load);
    for _ in 0..40 {
        sys.p_array = 0.5 * (lo + hi);
        if recharges(&sys, target) {
            hi = sys.p_array
        } else {
            lo = sys.p_array
        }
    }
    Some(hi)
}
