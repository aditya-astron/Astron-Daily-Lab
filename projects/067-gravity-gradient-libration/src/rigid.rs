// Astron Built this AP
// Project 067 | Attitude control | Gravity-gradient stability regions and eccentric-orbit pitch libration
// Author: Aditya (aditya-astron)
// Date: 2026-11-27 IST
//
// Full nonlinear rigid body in a Keplerian orbit under gravity-gradient torque.
// This is the reference model: the linear and planar equations in the other
// modules are only trusted once they agree with it.
//
// Units: mu = 1, a = 1, so the mean motion is 1 and one orbit is 2*pi.
// Orbit (LVLH) frame: x along-track, y = -orbit normal, z = nadir.
// Body axes are nominally aligned with it: 1 = roll, 2 = pitch, 3 = yaw.
//
// State (8): q0 q1 q2 q3 (C_bo, orbit -> body), w1 w2 w3 (inertial rate in
// body axes), nu (true anomaly).
//
//   C(q)  = (q0^2 - qv.qv) I + 2 qv qv^T - 2 q0 [qv x]
//   q0'   = -1/2 qv . w_bo,   qv' = 1/2 (q0 w_bo + qv x w_bo)
//   w_bo  = w - C (0, -nu', 0)            relative rate, body axes
//   I w'  = 3 mu/r^3 (c x I c) - w x I w, c = C (0,0,1) nadir in body
//   nu'   = sqrt(mu p) / r^2,  r = p / (1 + e cos nu)

pub type M3 = [[f64; 3]; 3];
pub type V3 = [f64; 3];

pub fn dcm(q: &[f64]) -> M3 {
    let (q0, q1, q2, q3) = (q[0], q[1], q[2], q[3]);
    let s = q0 * q0 - (q1 * q1 + q2 * q2 + q3 * q3);
    let v = [q1, q2, q3];
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = 2.0 * v[i] * v[j] + if i == j { s } else { 0.0 };
        }
    }
    // -2 q0 [qv x]
    c[0][1] += 2.0 * q0 * q3;
    c[0][2] -= 2.0 * q0 * q2;
    c[1][0] -= 2.0 * q0 * q3;
    c[1][2] += 2.0 * q0 * q1;
    c[2][0] += 2.0 * q0 * q2;
    c[2][1] -= 2.0 * q0 * q1;
    c
}

pub fn mv(m: &M3, v: &V3) -> V3 {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

pub fn cross(a: &V3, b: &V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Quaternion for C_bo = R1(roll) R2(pitch) R3(yaw), the 3-2-1 sequence.
pub fn quat_from_euler(roll: f64, pitch: f64, yaw: f64) -> [f64; 4] {
    let (cr, sr) = ((roll / 2.0).cos(), (roll / 2.0).sin());
    let (cp, sp) = ((pitch / 2.0).cos(), (pitch / 2.0).sin());
    let (cy, sy) = ((yaw / 2.0).cos(), (yaw / 2.0).sin());
    [
        cr * cp * cy + sr * sp * sy,
        sr * cp * cy - cr * sp * sy,
        cr * sp * cy + sr * cp * sy,
        cr * cp * sy - sr * sp * cy,
    ]
}

/// (roll, pitch, yaw) from C_bo for the same 3-2-1 sequence.
pub fn euler_from_dcm(c: &M3) -> (f64, f64, f64) {
    let pitch = (-c[0][2]).clamp(-1.0, 1.0).asin();
    let roll = c[1][2].atan2(c[2][2]);
    let yaw = c[0][1].atan2(c[0][0]);
    (roll, pitch, yaw)
}

#[derive(Clone, Copy, Debug)]
pub struct Body {
    pub i: V3, // principal inertias I1 (roll), I2 (pitch), I3 (yaw)
    pub e: f64,
}

impl Body {
    /// Normalised inertias from the Hughes ratios, I2 = 1:
    ///   k1 = (I2 - I3)/I1,  k3 = (I2 - I1)/I3
    pub fn from_k(k1: f64, k3: f64, e: f64) -> Body {
        let d = 1.0 - k1 * k3;
        Body {
            i: [(1.0 - k3) / d, 1.0, (1.0 - k1) / d],
            e,
        }
    }

    pub fn k1(&self) -> f64 {
        (self.i[1] - self.i[2]) / self.i[0]
    }
    pub fn k2(&self) -> f64 {
        (self.i[0] - self.i[2]) / self.i[1]
    }
    pub fn k3(&self) -> f64 {
        (self.i[1] - self.i[0]) / self.i[2]
    }

    fn p(&self) -> f64 {
        1.0 - self.e * self.e
    }

    pub fn nu_dot(&self, nu: f64) -> f64 {
        let p = self.p();
        let r = p / (1.0 + self.e * nu.cos());
        p.sqrt() / (r * r)
    }

    /// Initial state from Euler angles and their rates relative to the
    /// orbit frame. Angle rates are converted through the 3-2-1 kinematics,
    /// so a planar start (roll = yaw = 0) is exact for any pitch.
    pub fn state(&self, nu: f64, ang: (f64, f64, f64), rates: (f64, f64, f64)) -> Vec<f64> {
        let (r, p, _y) = ang;
        let (rd, pd, yd) = rates;
        let q = quat_from_euler(ang.0, ang.1, ang.2);
        let c = dcm(&q);
        // body-relative rate for C = R1(r) R2(p) R3(y)
        let w_bo = [
            rd - yd * p.sin(),
            pd * r.cos() + yd * r.sin() * p.cos(),
            -pd * r.sin() + yd * r.cos() * p.cos(),
        ];
        let w_oi = mv(&c, &[0.0, -self.nu_dot(nu), 0.0]);
        vec![
            q[0],
            q[1],
            q[2],
            q[3],
            w_bo[0] + w_oi[0],
            w_bo[1] + w_oi[1],
            w_bo[2] + w_oi[2],
            nu,
        ]
    }

    pub fn rhs(&self, _t: f64, y: &[f64], dy: &mut [f64]) {
        let c = dcm(&y[0..4]);
        let w = [y[4], y[5], y[6]];
        let nu = y[7];
        let nud = self.nu_dot(nu);
        let r = self.p() / (1.0 + self.e * nu.cos());

        let w_oi = mv(&c, &[0.0, -nud, 0.0]);
        let wr = [w[0] - w_oi[0], w[1] - w_oi[1], w[2] - w_oi[2]];
        let (q0, qv) = (y[0], [y[1], y[2], y[3]]);
        let qxw = cross(&qv, &wr);
        dy[0] = -0.5 * (qv[0] * wr[0] + qv[1] * wr[1] + qv[2] * wr[2]);
        for k in 0..3 {
            dy[1 + k] = 0.5 * (q0 * wr[k] + qxw[k]);
        }

        let nadir = [c[0][2], c[1][2], c[2][2]];
        let i = self.i;
        let ic = [i[0] * nadir[0], i[1] * nadir[1], i[2] * nadir[2]];
        let gg = cross(&nadir, &ic);
        let iw = [i[0] * w[0], i[1] * w[1], i[2] * w[2]];
        let gyro = cross(&w, &iw);
        let g = 3.0 / (r * r * r);
        for k in 0..3 {
            dy[4 + k] = (g * gg[k] - gyro[k]) / i[k];
        }
        dy[7] = nud;
    }

    /// Jacobi integral for a circular orbit (conserved only when e = 0):
    ///   H = 1/2 w_bo^T I w_bo + 3/2 c^T I c - 1/2 n^T I n
    /// with c the nadir and n the orbit normal in body axes. The potential
    /// part has a strict minimum at zero attitude only when I2 > I1 > I3,
    /// which is the Lagrange region.
    pub fn jacobi(&self, y: &[f64]) -> f64 {
        let c = dcm(&y[0..4]);
        let w_oi = mv(&c, &[0.0, -1.0, 0.0]);
        let wr = [y[4] - w_oi[0], y[5] - w_oi[1], y[6] - w_oi[2]];
        let nadir = [c[0][2], c[1][2], c[2][2]];
        let normal = [c[0][1], c[1][1], c[2][1]];
        let i = self.i;
        let quad = |v: &V3| i[0] * v[0] * v[0] + i[1] * v[1] * v[1] + i[2] * v[2] * v[2];
        0.5 * quad(&wr) + 1.5 * quad(&nadir) - 0.5 * quad(&normal)
    }

    pub fn angles(&self, y: &[f64]) -> (f64, f64, f64) {
        let n = (y[0] * y[0] + y[1] * y[1] + y[2] * y[2] + y[3] * y[3]).sqrt();
        let q = [y[0] / n, y[1] / n, y[2] / n, y[3] / n];
        euler_from_dcm(&dcm(&q))
    }
}
