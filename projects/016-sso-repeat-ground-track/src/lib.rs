// Astron Built this AP
// Project 016 | Orbital mechanics | Sun-synchronous repeat ground track design
// Author: Aditya (aditya-astron)
// Date: 2026-11-05 IST
//
// Design of circular sun-synchronous orbits whose ground track repeats exactly
// after N revolutions in D days, using first-order J2 secular theory, plus a
// numerical J2 propagator to check the analytic design.

use std::f64::consts::PI;

pub const MU: f64 = 398_600.441_8; // km^3/s^2, WGS-84
pub const RE: f64 = 6_378.137; // km, equatorial radius
pub const J2: f64 = 1.082_626_68e-3;
pub const OMEGA_E: f64 = 7.2921158553e-5; // rad/s, sidereal rotation
pub const TROPICAL_YEAR_S: f64 = 365.242_189_7 * 86_400.0;
/// Mean motion of the mean Sun: the required nodal precession for SSO.
pub const SUN_RATE: f64 = 2.0 * PI / TROPICAL_YEAR_S;

/// First-order J2 secular rates for a circular orbit (e = 0), in rad/s.
#[derive(Debug, Clone, Copy)]
pub struct SecularRates {
    pub raan_dot: f64,
    pub argp_dot: f64,
    pub mean_anom_dot: f64,
}

impl SecularRates {
    /// Rate of the argument of latitude u = w + M; one nodal period is 2*pi / u_dot.
    pub fn arg_lat_dot(&self) -> f64 {
        self.argp_dot + self.mean_anom_dot
    }
}

/// Kozai/Brouwer first-order secular rates with e = 0:
///   dRAAN/dt = -3/2 n J2 (Re/a)^2 cos i
///   dw/dt    =  3/4 n J2 (Re/a)^2 (5 cos^2 i - 1)
///   dM/dt    =  n [1 + 3/4 J2 (Re/a)^2 (3 cos^2 i - 1)]
pub fn secular_rates(a: f64, inc: f64) -> SecularRates {
    let n = (MU / a.powi(3)).sqrt();
    let k = J2 * (RE / a).powi(2);
    let c2 = inc.cos().powi(2);
    SecularRates {
        raan_dot: -1.5 * n * k * inc.cos(),
        argp_dot: 0.75 * n * k * (5.0 * c2 - 1.0),
        mean_anom_dot: n * (1.0 + 0.75 * k * (3.0 * c2 - 1.0)),
    }
}

/// Inclination that makes the node precess with the mean Sun at semi-major axis a.
/// Returns None above the SSO altitude limit where |cos i| would exceed 1.
pub fn sso_inclination(a: f64) -> Option<f64> {
    let n = (MU / a.powi(3)).sqrt();
    let k = J2 * (RE / a).powi(2);
    let ci = -SUN_RATE / (1.5 * n * k);
    if ci.abs() <= 1.0 {
        Some(ci.acos())
    } else {
        None
    }
}

/// Length of a nodal day: time for Earth to turn once relative to the orbit plane.
/// For an SSO this equals the mean solar day to within the tropical/sidereal constants.
pub fn nodal_day(raan_dot: f64) -> f64 {
    2.0 * PI / (OMEGA_E - raan_dot)
}

#[derive(Debug, Clone, Copy)]
pub struct RepeatOrbit {
    pub revs: u32,
    pub days: u32,
    pub a: f64,
    pub inc: f64,
    pub nodal_period: f64,
    pub nodal_day: f64,
}

impl RepeatOrbit {
    pub fn altitude(&self) -> f64 {
        self.a - RE
    }

    /// Spacing between adjacent ascending-node tracks at the equator once the full
    /// cycle has been flown (km).
    pub fn equatorial_spacing(&self) -> f64 {
        2.0 * PI * RE / self.revs as f64
    }

    /// Westward shift of the node between consecutive revolutions (km at the equator).
    pub fn rev_to_rev_shift(&self) -> f64 {
        2.0 * PI * RE * self.days as f64 / self.revs as f64
    }

    /// Days between passes over adjacent tracks: the smallest d with
    /// d*N = +/-1 (mod D). Returns (d, true) if the adjacent track lies east
    /// of the reference track, (d, false) if west.
    pub fn subcycle(&self) -> (u32, bool) {
        subcycle(self.revs, self.days)
    }
}

pub fn subcycle(revs: u32, days: u32) -> (u32, bool) {
    if days == 1 {
        return (1, true);
    }
    let r = revs % days;
    for d in 1..=days {
        let m = (d * r) % days;
        // Rev k has node longitude -(2*pi/N)(k*D mod N). Rev k = (d*N + s)/D
        // occurs s/N day after d whole days and sits s track spacings west of
        // the reference, so d*N = -s (mod D) with s = +1 (west) or -1 (east).
        // For D = 2 both hold; the east neighbour comes first (d - 1/N days).
        if m == 1 {
            return (d, true);
        }
        if m == days - 1 {
            return (d, false);
        }
    }
    (days, true)
}

/// Residual of the repeat condition at semi-major axis a with the SSO inclination:
///   N * T_nodal - D * T_nodal_day = 0, expressed as a rate difference.
fn repeat_residual(a: f64, revs: u32, days: u32) -> Option<f64> {
    let inc = sso_inclination(a)?;
    let r = secular_rates(a, inc);
    let q = revs as f64 / days as f64;
    Some(r.arg_lat_dot() - q * (OMEGA_E - r.raan_dot))
}

/// Solve for the SSO with N revs in D nodal days. The residual is monotonic
/// decreasing in a over the LEO range, so a bracketed bisection is robust.
pub fn solve_repeat_sso(revs: u32, days: u32) -> Option<RepeatOrbit> {
    let mut lo = RE + 150.0;
    let mut hi = RE + 5_500.0;
    let flo = repeat_residual(lo, revs, days)?;
    let fhi = repeat_residual(hi, revs, days)?;
    if flo.signum() == fhi.signum() {
        return None;
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        let fm = repeat_residual(mid, revs, days)?;
        if fm.signum() == flo.signum() {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo < 1e-10 {
            break;
        }
    }
    let a = 0.5 * (lo + hi);
    let inc = sso_inclination(a)?;
    let r = secular_rates(a, inc);
    Some(RepeatOrbit {
        revs,
        days,
        a,
        inc,
        nodal_period: 2.0 * PI / r.arg_lat_dot(),
        nodal_day: nodal_day(r.raan_dot),
    })
}

fn gcd(mut a: u32, mut b: u32) -> u32 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

/// All coprime (N, D) repeat SSOs with D <= max_days inside an altitude band.
/// N/D is bracketed with the Keplerian mean motion at the band edges, padded
/// slightly to cover the J2 correction.
pub fn enumerate_repeat_orbits(h_min: f64, h_max: f64, max_days: u32) -> Vec<RepeatOrbit> {
    let revs_per_day = |h: f64| {
        let a: f64 = RE + h;
        (MU / a.powi(3)).sqrt() / (OMEGA_E - SUN_RATE)
    };
    let q_hi = revs_per_day(h_min - 30.0);
    let q_lo = revs_per_day(h_max + 30.0);
    let mut out = Vec::new();
    for d in 1..=max_days {
        let n_lo = (q_lo * d as f64).floor() as u32;
        let n_hi = (q_hi * d as f64).ceil() as u32;
        for n in n_lo..=n_hi {
            if gcd(n, d) != 1 {
                continue;
            }
            if let Some(o) = solve_repeat_sso(n, d) {
                let h = o.altitude();
                if h >= h_min && h <= h_max {
                    out.push(o);
                }
            }
        }
    }
    out.sort_by(|x, y| x.a.partial_cmp(&y.a).unwrap());
    out
}

// ---------------------------------------------------------------------------
// Numerical check: Cartesian two-body + J2 propagation (ECI, RK4).

pub type State = [f64; 6];

fn deriv(s: &State) -> State {
    let (x, y, z) = (s[0], s[1], s[2]);
    let r2 = x * x + y * y + z * z;
    let r = r2.sqrt();
    let z2 = z * z / r2;
    let kep = -MU / (r2 * r);
    // J2 acceleration: 3/2 J2 mu Re^2 / r^5 * [x(5z^2/r^2-1), y(5z^2/r^2-1), z(5z^2/r^2-3)]
    let c = 1.5 * J2 * MU * RE * RE / (r2 * r2 * r);
    [
        s[3],
        s[4],
        s[5],
        kep * x + c * x * (5.0 * z2 - 1.0),
        kep * y + c * y * (5.0 * z2 - 1.0),
        kep * z + c * z * (5.0 * z2 - 3.0),
    ]
}

pub fn rk4_step(s: &State, h: f64) -> State {
    let add = |a: &State, b: &State, k: f64| {
        let mut o = [0.0; 6];
        for i in 0..6 {
            o[i] = a[i] + k * b[i];
        }
        o
    };
    let k1 = deriv(s);
    let k2 = deriv(&add(s, &k1, 0.5 * h));
    let k3 = deriv(&add(s, &k2, 0.5 * h));
    let k4 = deriv(&add(s, &k3, h));
    let mut o = [0.0; 6];
    for i in 0..6 {
        o[i] = s[i] + h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
    o
}

/// Circular orbit started at the ascending node on the inertial x-axis.
pub fn state_at_node(a: f64, inc: f64) -> State {
    let v = (MU / a).sqrt();
    [a, 0.0, 0.0, 0.0, v * inc.cos(), v * inc.sin()]
}

/// Osculating semi-major axis and inclination.
pub fn osc_a_inc(s: &State) -> (f64, f64) {
    let r = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt();
    let v2 = s[3] * s[3] + s[4] * s[4] + s[5] * s[5];
    let a = 1.0 / (2.0 / r - v2 / MU);
    let hx = s[1] * s[5] - s[2] * s[4];
    let hy = s[2] * s[3] - s[0] * s[5];
    let hz = s[0] * s[4] - s[1] * s[3];
    let h = (hx * hx + hy * hy + hz * hz).sqrt();
    (a, (hz / h).acos())
}

/// Orbit-averaged osculating a and i over one nodal period, sampled uniformly
/// (a periodic integrand, so the rectangle rule converges very fast).
pub fn orbit_averaged(a0: f64, i0: f64, period: f64, steps: usize) -> (f64, f64) {
    let h = period / steps as f64;
    let mut s = state_at_node(a0, i0);
    let (mut sa, mut si) = (0.0, 0.0);
    for _ in 0..steps {
        let (a, i) = osc_a_inc(&s);
        sa += a;
        si += i;
        s = rk4_step(&s, h);
    }
    (sa / steps as f64, si / steps as f64)
}

/// Find the osculating (a, i) at the ascending node whose orbit average equals
/// the requested mean (a, i). Short-period J2 terms put the osculating a about
/// 9 km above the mean value at the node of a near-polar LEO.
pub fn mean_to_osculating_at_node(a_mean: f64, i_mean: f64, period: f64) -> (f64, f64) {
    let (mut a0, mut i0) = (a_mean, i_mean);
    for _ in 0..6 {
        let (am, im) = orbit_averaged(a0, i0, period, 4000);
        a0 += a_mean - am;
        i0 += i_mean - im;
    }
    (a0, i0)
}

#[derive(Debug, Clone)]
pub struct NodeCrossing {
    pub t: f64,
    /// Right ascension of the node in the inertial frame (rad, unwrapped later).
    pub ra: f64,
}

/// Propagate and record ascending node crossings (z goes - to +). Each crossing
/// is located by bisection on a re-integrated sub-step.
pub fn propagate_nodes(s0: State, duration: f64, h: f64) -> Vec<NodeCrossing> {
    let mut out = Vec::new();
    let mut s = s0;
    let mut t = 0.0;
    let steps = (duration / h).ceil() as usize;
    for _ in 0..steps {
        let next = rk4_step(&s, h);
        if s[2] < 0.0 && next[2] >= 0.0 {
            let (mut lo, mut hi) = (0.0, h);
            let mut sc = next;
            for _ in 0..60 {
                let mid = 0.5 * (lo + hi);
                sc = rk4_step(&s, mid);
                if sc[2] < 0.0 {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            out.push(NodeCrossing {
                t: t + 0.5 * (lo + hi),
                ra: sc[1].atan2(sc[0]),
            });
        }
        s = next;
        t += h;
    }
    out
}

/// Least-squares slope of y against x.
pub fn slope(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let mut sxy = 0.0;
    let mut sxx = 0.0;
    for (xi, yi) in x.iter().zip(y) {
        sxy += (xi - mx) * (yi - my);
        sxx += (xi - mx) * (xi - mx);
    }
    sxy / sxx
}

pub fn unwrap(angles: &[f64]) -> Vec<f64> {
    let mut out = Vec::with_capacity(angles.len());
    let mut off = 0.0;
    for (k, &a) in angles.iter().enumerate() {
        if k > 0 {
            let d = a + off - out[k - 1];
            if d > PI {
                off -= 2.0 * PI;
            } else if d < -PI {
                off += 2.0 * PI;
            }
        }
        out.push(a + off);
    }
    out
}

/// Wrap an angle to (-pi, pi].
pub fn wrap_pi(x: f64) -> f64 {
    let mut y = x % (2.0 * PI);
    if y > PI {
        y -= 2.0 * PI;
    } else if y <= -PI {
        y += 2.0 * PI;
    }
    y
}

#[derive(Debug, Clone, Copy)]
pub struct PropagationCheck {
    pub nodal_period: f64,
    pub raan_rate: f64,
    /// Earth-fixed longitude of the ascending node after N revolutions,
    /// relative to the first node (km along the equator).
    pub closure_km: f64,
}

/// Fly a full repeat cycle numerically and measure what the design promised.
pub fn verify_by_propagation(o: &RepeatOrbit, h: f64) -> PropagationCheck {
    let (a0, i0) = mean_to_osculating_at_node(o.a, o.inc, o.nodal_period);
    let s0 = state_at_node(a0, i0);
    let duration = o.revs as f64 * o.nodal_period + 0.5 * o.nodal_period;
    let mut nodes = vec![NodeCrossing { t: 0.0, ra: 0.0 }];
    nodes.extend(propagate_nodes(s0, duration, h));
    let t: Vec<f64> = nodes.iter().map(|n| n.t).collect();
    let ra = unwrap(&nodes.iter().map(|n| n.ra).collect::<Vec<_>>());
    let revs = o.revs as usize;
    assert!(nodes.len() > revs, "propagation ended before the cycle closed");
    // Earth-fixed node longitude: inertial RA minus Earth rotation (GMST = 0 at t = 0).
    let lon_end = wrap_pi(ra[revs] - OMEGA_E * t[revs]);
    PropagationCheck {
        nodal_period: (t[revs] - t[0]) / o.revs as f64,
        raan_rate: slope(&t, &ra),
        closure_km: lon_end * RE,
    }
}
