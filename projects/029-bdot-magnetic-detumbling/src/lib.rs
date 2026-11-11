// Astron Built this AP
// Project 029 | Attitude control | B-dot magnetic detumbling of a 3U CubeSat
// Author: Aditya (aditya-astron)
// Date: 2026-11-11 IST

//! Rigid-body attitude simulation of a magnetorquer-only detumble.
//!
//! Frames: ECI (inertial), ECEF (Earth-fixed, for the dipole axis), body.
//! The quaternion q = [w, x, y, z] maps body vectors into ECI.

pub type V3 = [f64; 3];
pub type M3 = [[f64; 3]; 3];

pub const MU_EARTH: f64 = 3.986004418e14; // m^3/s^2
pub const R_EARTH: f64 = 6.3712e6; // IGRF reference radius, m
pub const OMEGA_EARTH: f64 = 7.2921159e-5; // rad/s

// ---------------------------------------------------------------- vectors

pub fn add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
pub fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
pub fn scale(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
pub fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
pub fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
pub fn norm(a: V3) -> f64 {
    dot(a, a).sqrt()
}
pub fn unit(a: V3) -> V3 {
    scale(a, 1.0 / norm(a))
}

// ------------------------------------------------------------ quaternions

pub type Quat = [f64; 4];

pub fn qmul(p: Quat, q: Quat) -> Quat {
    [
        p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3],
        p[0] * q[1] + p[1] * q[0] + p[2] * q[3] - p[3] * q[2],
        p[0] * q[2] - p[1] * q[3] + p[2] * q[0] + p[3] * q[1],
        p[0] * q[3] + p[1] * q[2] - p[2] * q[1] + p[3] * q[0],
    ]
}

/// Body -> ECI.
pub fn rotate(q: Quat, v: V3) -> V3 {
    let r = qmul(qmul(q, [0.0, v[0], v[1], v[2]]), [q[0], -q[1], -q[2], -q[3]]);
    [r[1], r[2], r[3]]
}

/// ECI -> body.
pub fn rotate_inv(q: Quat, v: V3) -> V3 {
    rotate([q[0], -q[1], -q[2], -q[3]], v)
}

fn qnormalize(q: Quat) -> Quat {
    let n = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
    [q[0] / n, q[1] / n, q[2] / n, q[3] / n]
}

// ------------------------------------------------------- orbit and field

/// Circular Keplerian orbit; J2 and drag do not matter over a few orbits of
/// detumbling.
#[derive(Clone, Copy, Debug)]
pub struct Orbit {
    pub radius: f64,
    pub incl: f64,
    pub raan: f64,
    pub u0: f64,
}

impl Orbit {
    pub fn mean_motion(&self) -> f64 {
        (MU_EARTH / self.radius.powi(3)).sqrt()
    }
    pub fn period(&self) -> f64 {
        2.0 * std::f64::consts::PI / self.mean_motion()
    }
    pub fn position(&self, t: f64) -> V3 {
        let u = self.u0 + self.mean_motion() * t;
        let (su, cu) = u.sin_cos();
        let (si, ci) = self.incl.sin_cos();
        let (so, co) = self.raan.sin_cos();
        scale(
            [co * cu - so * su * ci, so * cu + co * su * ci, su * si],
            self.radius,
        )
    }
}

/// Centred dipole: B = B0 (Re/r)^3 [3 (m.r) r - m], m the unit dipole axis.
/// With IGRF degree-1 coefficients the axis in ECEF is (g11, h11, g10)/B0,
/// which points roughly at the geographic south pole.
#[derive(Clone, Copy, Debug)]
pub struct Dipole {
    pub b0: f64,
    pub axis_ecef: V3,
    pub rotating: bool,
    pub theta0: f64,
}

impl Dipole {
    /// IGRF-13 epoch 2020 degree-1 coefficients (nT).
    pub fn igrf2020() -> Dipole {
        let g = [-1450.9e-9, 4652.5e-9, -29404.8e-9];
        let b0 = norm(g);
        Dipole { b0, axis_ecef: scale(g, 1.0 / b0), rotating: true, theta0: 0.0 }
    }

    /// Axis on -z, Earth not rotating: the setting of the averaging theory.
    pub fn aligned(b0: f64) -> Dipole {
        Dipole { b0, axis_ecef: [0.0, 0.0, -1.0], rotating: false, theta0: 0.0 }
    }

    pub fn axis_eci(&self, t: f64) -> V3 {
        let th = self.theta0 + if self.rotating { OMEGA_EARTH * t } else { 0.0 };
        let (s, c) = th.sin_cos();
        let m = self.axis_ecef;
        [c * m[0] - s * m[1], s * m[0] + c * m[1], m[2]]
    }

    pub fn field(&self, r: V3, t: f64) -> V3 {
        let rn = norm(r);
        let rh = scale(r, 1.0 / rn);
        let m = self.axis_eci(t);
        let k = self.b0 * (R_EARTH / rn).powi(3);
        scale(sub(scale(rh, 3.0 * dot(m, rh)), m), k)
    }
}

// --------------------------------------------------------- spacecraft

#[derive(Clone, Copy, Debug)]
pub struct Spacecraft {
    /// Principal moments of inertia, kg m^2.
    pub inertia: V3,
    /// Per-axis magnetorquer saturation, A m^2.
    pub m_max: f64,
}

impl Spacecraft {
    pub fn kinetic_energy(&self, w: V3) -> f64 {
        0.5 * (self.inertia[0] * w[0] * w[0]
            + self.inertia[1] * w[1] * w[1]
            + self.inertia[2] * w[2] * w[2])
    }
    pub fn j_min(&self) -> f64 {
        self.inertia.iter().cloned().fold(f64::INFINITY, f64::min)
    }
}

/// Control law.
///
/// `Ideal`: m = (k/|B|) w x b, so the torque is -k (I - b b^T) w with no
/// saturation and perfect rate knowledge. This is the law the averaging
/// analysis is about.
///
/// `Sampled`: what flight software does. Each `period` the torquers are
/// switched off for `quiet` seconds around a magnetometer read, B-dot is the
/// backward difference of two reads, and m = -(k/|B|^2) Bdot is held, scaled
/// down as a vector if any axis exceeds m_max. Since Bdot_body ~ -w x B, the
/// two laws agree when the spin is fast compared with the orbit rate.
#[derive(Clone, Copy, Debug)]
pub enum Law {
    Ideal { k: f64 },
    Sampled { k: f64, period: f64, quiet: f64, noise: f64 },
}

#[derive(Clone, Copy, Debug)]
pub struct Scenario {
    pub sc: Spacecraft,
    pub orbit: Orbit,
    pub field: Dipole,
    pub law: Law,
    pub w0: V3,
    pub q0: Quat,
    pub dt: f64,
    pub seed: u64,
}

#[derive(Clone, Copy, Debug)]
pub struct Sample {
    pub t: f64,
    pub w: V3,
    pub energy: f64,
    /// Sampled law: largest axis of the held dipole command (A m^2).
    /// Ideal law: magnitude of the dipole it would need.
    pub dipole: f64,
}

/// Avanzini & Giulietti (JGCD 2012) gain for the normalised law:
/// k* = 2 n (1 + sin xi) J_min, xi the orbit inclination to the geomagnetic
/// equator. Above this the decay rate is limited by how fast the field
/// direction sweeps round, not by the gain.
pub fn avanzini_gain(n: f64, xi: f64, j_min: f64) -> f64 {
    2.0 * n * (1.0 + xi.sin()) * j_min
}

/// Saturate as a vector so the dipole direction is preserved.
pub fn saturate(m: V3, m_max: f64) -> V3 {
    let peak = m.iter().fold(0.0f64, |a, x| a.max(x.abs()));
    if peak > m_max { scale(m, m_max / peak) } else { m }
}

// xorshift64* plus Box-Muller: reproducible sensor noise without a crate.
struct Rng(u64);
impl Rng {
    fn uniform(&mut self) -> f64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        let x = self.0.wrapping_mul(0x2545F4914F6CDD1D);
        ((x >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    }
    fn gauss(&mut self) -> f64 {
        let (u1, u2) = (self.uniform(), self.uniform());
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

#[derive(Clone, Copy)]
struct State {
    q: Quat,
    w: V3,
}

/// Body-frame torque from a body-frame dipole, or from the ideal law.
fn torque(s: &Scenario, x: &State, t: f64, m_body: V3) -> V3 {
    let b_eci = s.field.field(s.orbit.position(t), t);
    let b = rotate_inv(x.q, b_eci);
    match s.law {
        Law::Ideal { k } => {
            let bh = unit(b);
            let m = scale(cross(x.w, bh), k / norm(b));
            cross(m, b)
        }
        Law::Sampled { .. } => cross(m_body, b),
    }
}

fn deriv(s: &Scenario, x: &State, t: f64, m_body: V3) -> State {
    let j = s.sc.inertia;
    let tau = torque(s, x, t, m_body);
    let h = [j[0] * x.w[0], j[1] * x.w[1], j[2] * x.w[2]];
    let gyro = cross(x.w, h);
    let wd = [
        (tau[0] - gyro[0]) / j[0],
        (tau[1] - gyro[1]) / j[1],
        (tau[2] - gyro[2]) / j[2],
    ];
    let qd = qmul(x.q, [0.0, x.w[0], x.w[1], x.w[2]]);
    State { q: [0.5 * qd[0], 0.5 * qd[1], 0.5 * qd[2], 0.5 * qd[3]], w: wd }
}

fn axpy(x: &State, d: &State, h: f64) -> State {
    State {
        q: [x.q[0] + h * d.q[0], x.q[1] + h * d.q[1], x.q[2] + h * d.q[2], x.q[3] + h * d.q[3]],
        w: add(x.w, scale(d.w, h)),
    }
}

fn rk4(s: &Scenario, x: &State, t: f64, h: f64, m: V3) -> State {
    let k1 = deriv(s, x, t, m);
    let k2 = deriv(s, &axpy(x, &k1, 0.5 * h), t + 0.5 * h, m);
    let k3 = deriv(s, &axpy(x, &k2, 0.5 * h), t + 0.5 * h, m);
    let k4 = deriv(s, &axpy(x, &k3, h), t + h, m);
    let mut q = [0.0; 4];
    for i in 0..4 {
        q[i] = x.q[i] + h / 6.0 * (k1.q[i] + 2.0 * k2.q[i] + 2.0 * k3.q[i] + k4.q[i]);
    }
    let mut w = [0.0; 3];
    for i in 0..3 {
        w[i] = x.w[i] + h / 6.0 * (k1.w[i] + 2.0 * k2.w[i] + 2.0 * k3.w[i] + k4.w[i]);
    }
    State { q: qnormalize(q), w }
}

/// Integrate for `t_end` seconds, logging every `log_every` seconds.
/// For the sampled law, `period` and `quiet` must be multiples of `dt`.
pub fn simulate(s: &Scenario, t_end: f64, log_every: f64) -> Vec<Sample> {
    let mut x = State { q: s.q0, w: s.w0 };
    let mut rng = Rng(s.seed.max(1));
    let steps = (t_end / s.dt).round() as usize;
    let log_stride = ((log_every / s.dt).round() as usize).max(1);

    let (per_steps, quiet_steps) = match s.law {
        Law::Sampled { period, quiet, .. } => (
            (period / s.dt).round() as usize,
            (quiet / s.dt).round() as usize,
        ),
        Law::Ideal { .. } => (0, 0),
    };

    let mut b_prev: Option<V3> = None;
    let mut m_hold = [0.0; 3];
    let mut out = Vec::with_capacity(steps / log_stride + 2);

    for i in 0..=steps {
        let t = i as f64 * s.dt;
        let mut m_now = [0.0; 3];

        if let Law::Sampled { k, period, noise, .. } = s.law {
            let phase = i % per_steps;
            if phase == 0 {
                // Torquers have been off for `quiet` s: read the magnetometer.
                let b_true = rotate_inv(x.q, s.field.field(s.orbit.position(t), t));
                let b_meas = [
                    b_true[0] + noise * rng.gauss(),
                    b_true[1] + noise * rng.gauss(),
                    b_true[2] + noise * rng.gauss(),
                ];
                m_hold = match b_prev {
                    Some(bp) => {
                        let bdot = scale(sub(b_meas, bp), 1.0 / period);
                        let b2 = dot(b_meas, b_meas);
                        saturate(scale(bdot, -k / b2), s.sc.m_max)
                    }
                    None => [0.0; 3],
                };
                b_prev = Some(b_meas);
            }
            // On between the quiet windows that bracket each read.
            if phase >= quiet_steps && phase < per_steps - quiet_steps {
                m_now = m_hold;
            }
        }

        if i % log_stride == 0 {
            let dip = match s.law {
                Law::Ideal { k } => {
                    let b = rotate_inv(x.q, s.field.field(s.orbit.position(t), t));
                    k * norm(cross(x.w, unit(b))) / norm(b)
                }
                // Largest commanded axis, whether or not this instant is quiet.
                Law::Sampled { .. } => m_hold.iter().fold(0.0f64, |a, x| a.max(x.abs())),
            };
            out.push(Sample { t, w: x.w, energy: s.sc.kinetic_energy(x.w), dipole: dip });
        }
        if i == steps {
            break;
        }
        x = rk4(s, &x, t, s.dt, m_now);
    }
    out
}

/// First time after which |w| stays below `thresh` for the rest of the run.
pub fn settle_time(log: &[Sample], thresh: f64) -> Option<f64> {
    let mut idx = None;
    for (i, smp) in log.iter().enumerate().rev() {
        if norm(smp.w) >= thresh {
            break;
        }
        idx = Some(i);
    }
    idx.map(|i| log[i].t)
}

// ----------------------------------------------------- averaging theory

/// Orbit average of b b^T for a circular orbit inclined `xi` to the
/// equator of an aligned, non-rotating dipole (orbit-plane coordinates:
/// x to the node, z along the dipole axis). Midpoint rule on a periodic
/// integrand converges geometrically.
pub fn averaged_bbt(xi: f64, n_quad: usize) -> M3 {
    let orbit = Orbit { radius: R_EARTH + 500e3, incl: xi, raan: 0.0, u0: 0.0 };
    let dip = Dipole::aligned(3e-5);
    let mut m = [[0.0; 3]; 3];
    for k in 0..n_quad {
        let t = (k as f64 + 0.5) / n_quad as f64 * orbit.period();
        let b = unit(dip.field(orbit.position(t), t));
        for i in 0..3 {
            for j in 0..3 {
                m[i][j] += b[i] * b[j] / n_quad as f64;
            }
        }
    }
    m
}

/// Eigenvalues of a symmetric 3x3 matrix by cyclic Jacobi, ascending.
pub fn sym_eigenvalues(a: M3) -> V3 {
    let mut a = a;
    for _ in 0..50 {
        let off = a[0][1].powi(2) + a[0][2].powi(2) + a[1][2].powi(2);
        if off < 1e-30 {
            break;
        }
        for (p, q) in [(0, 1), (0, 2), (1, 2)] {
            if a[p][q].abs() < 1e-300 {
                continue;
            }
            let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
            let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
            let t = if theta == 0.0 { 1.0 } else { t };
            let c = 1.0 / (t * t + 1.0).sqrt();
            let s = t * c;
            // A <- J^T A J with the (p, q) rotation
            for k in 0..3 {
                let akp = a[k][p];
                let akq = a[k][q];
                a[k][p] = c * akp - s * akq;
                a[k][q] = s * akp + c * akq;
            }
            for k in 0..3 {
                let apk = a[p][k];
                let aqk = a[q][k];
                a[p][k] = c * apk - s * aqk;
                a[q][k] = s * apk + c * aqk;
            }
        }
    }
    let mut ev = [a[0][0], a[1][1], a[2][2]];
    ev.sort_by(|x, y| x.partial_cmp(y).unwrap());
    ev
}

/// Averaged model J dw/dt = -k (I - <bb^T>) w (isotropic J): the slowest
/// decay rate of |w| is (k/J)(1 - lambda_max).
pub fn averaged_slowest_rate(k: f64, j: f64, xi: f64) -> f64 {
    let ev = sym_eigenvalues(averaged_bbt(xi, 720));
    k / j * (1.0 - ev[2])
}
