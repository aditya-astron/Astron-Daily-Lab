// Astron Built this AP
// Project 005 | Attitude determination | Wahba's problem: TRIAD, Davenport q-method and QUEST
// Author: Aditya (aditya-astron)
// Date: 2026-11-01 IST
//
// Static attitude determination from vector observations.
//
// Conventions (Shuster / Markley & Crassidis):
//   * Attitude matrix A maps reference-frame vectors into the body frame: b = A r.
//   * Quaternion is scalar-last, q = [q1 q2 q3 q4], with
//       A(q) = (q4^2 - |qv|^2) I + 2 qv qv^T - 2 q4 [qv x]
//   * Small attitude errors are expressed in the body frame:
//       A_est A_true^T = I - [dtheta x]
//
// Wahba's problem: find the rotation A minimising
//   L(A) = 1/2 sum_i a_i |b_i - A r_i|^2
// which is the same as maximising tr(A B^T) with B = sum_i a_i b_i r_i^T.

pub type V3 = [f64; 3];
pub type M3 = [[f64; 3]; 3];
pub type Quat = [f64; 4];

// ---------- small 3-vector / 3x3 helpers (no external crates) ----------

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

pub fn scale(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

pub fn add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

pub fn normalize(a: V3) -> V3 {
    scale(a, 1.0 / norm(a))
}

pub fn outer(a: V3, b: V3) -> M3 {
    let mut m = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            m[i][j] = a[i] * b[j];
        }
    }
    m
}

pub fn mat_add(a: &M3, b: &M3) -> M3 {
    let mut m = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            m[i][j] = a[i][j] + b[i][j];
        }
    }
    m
}

pub fn mat_scale(a: &M3, s: f64) -> M3 {
    let mut m = *a;
    for row in m.iter_mut() {
        for x in row.iter_mut() {
            *x *= s;
        }
    }
    m
}

pub fn mat_mul(a: &M3, b: &M3) -> M3 {
    let mut m = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            m[i][j] = (0..3).map(|k| a[i][k] * b[k][j]).sum();
        }
    }
    m
}

pub fn transpose(a: &M3) -> M3 {
    let mut m = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            m[i][j] = a[j][i];
        }
    }
    m
}

pub fn mat_vec(a: &M3, v: V3) -> V3 {
    [dot(a[0], v), dot(a[1], v), dot(a[2], v)]
}

pub fn identity() -> M3 {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}

pub fn trace(a: &M3) -> f64 {
    a[0][0] + a[1][1] + a[2][2]
}

pub fn det3(a: &M3) -> f64 {
    a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
}

/// Adjugate (transpose of the cofactor matrix), so a * adj(a) = det(a) I.
pub fn adj3(a: &M3) -> M3 {
    let c = |r0: usize, r1: usize, c0: usize, c1: usize| a[r0][c0] * a[r1][c1] - a[r0][c1] * a[r1][c0];
    [
        [c(1, 2, 1, 2), -c(0, 2, 1, 2), c(0, 1, 1, 2)],
        [-c(1, 2, 0, 2), c(0, 2, 0, 2), -c(0, 1, 0, 2)],
        [c(1, 2, 0, 1), -c(0, 2, 0, 1), c(0, 1, 0, 1)],
    ]
}

pub fn inv3(a: &M3) -> M3 {
    mat_scale(&adj3(a), 1.0 / det3(a))
}

// ---------- quaternion / attitude matrix ----------

pub fn quat_normalize(q: Quat) -> Quat {
    let n = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
    [q[0] / n, q[1] / n, q[2] / n, q[3] / n]
}

pub fn quat_to_dcm(q: Quat) -> M3 {
    let (q1, q2, q3, q4) = (q[0], q[1], q[2], q[3]);
    let s = q4 * q4 - (q1 * q1 + q2 * q2 + q3 * q3);
    [
        [s + 2.0 * q1 * q1, 2.0 * (q1 * q2 + q4 * q3), 2.0 * (q1 * q3 - q4 * q2)],
        [2.0 * (q1 * q2 - q4 * q3), s + 2.0 * q2 * q2, 2.0 * (q2 * q3 + q4 * q1)],
        [2.0 * (q1 * q3 + q4 * q2), 2.0 * (q2 * q3 - q4 * q1), s + 2.0 * q3 * q3],
    ]
}

/// Shepperd's method: pick the largest of (q4, q1, q2, q3) to divide by, so
/// the extraction never divides by a near-zero component.
pub fn dcm_to_quat(a: &M3) -> Quat {
    let tr = trace(a);
    let cands = [a[0][0], a[1][1], a[2][2], tr];
    let k = (0..4).max_by(|&i, &j| cands[i].partial_cmp(&cands[j]).unwrap()).unwrap();
    let q = match k {
        3 => {
            let q4 = 0.5 * (1.0 + tr).sqrt();
            let f = 0.25 / q4;
            [f * (a[1][2] - a[2][1]), f * (a[2][0] - a[0][2]), f * (a[0][1] - a[1][0]), q4]
        }
        0 => {
            let q1 = 0.5 * (1.0 + 2.0 * a[0][0] - tr).sqrt();
            let f = 0.25 / q1;
            [q1, f * (a[0][1] + a[1][0]), f * (a[0][2] + a[2][0]), f * (a[1][2] - a[2][1])]
        }
        1 => {
            let q2 = 0.5 * (1.0 + 2.0 * a[1][1] - tr).sqrt();
            let f = 0.25 / q2;
            [f * (a[0][1] + a[1][0]), q2, f * (a[1][2] + a[2][1]), f * (a[2][0] - a[0][2])]
        }
        _ => {
            let q3 = 0.5 * (1.0 + 2.0 * a[2][2] - tr).sqrt();
            let f = 0.25 / q3;
            [f * (a[0][2] + a[2][0]), f * (a[1][2] + a[2][1]), q3, f * (a[0][1] - a[1][0])]
        }
    };
    quat_normalize(q)
}

/// Rotation of `angle` rad about unit `axis`, as an attitude quaternion.
pub fn quat_from_axis_angle(axis: V3, angle: f64) -> Quat {
    let e = normalize(axis);
    let s = (0.5 * angle).sin();
    [e[0] * s, e[1] * s, e[2] * s, (0.5 * angle).cos()]
}

/// Small-angle error vector in the body frame: A_est A_true^T = I - [dtheta x].
pub fn error_vector(a_est: &M3, a_true: &M3) -> V3 {
    let d = mat_mul(a_est, &transpose(a_true));
    [
        0.5 * (d[1][2] - d[2][1]),
        0.5 * (d[2][0] - d[0][2]),
        0.5 * (d[0][1] - d[1][0]),
    ]
}

/// Exact rotation angle between two attitudes (rad), valid for any size.
pub fn error_angle(a_est: &M3, a_true: &M3) -> f64 {
    // atan2 of (sin, cos) parts keeps full precision near 0 and near pi,
    // where acos of the trace alone bottoms out around 1e-8 rad.
    let d = mat_mul(a_est, &transpose(a_true));
    let s = norm([d[1][2] - d[2][1], d[2][0] - d[0][2], d[0][1] - d[1][0]]) * 0.5;
    s.atan2((trace(&d) - 1.0) * 0.5)
}

// ---------- the problem ----------

/// One vector observation: body-frame measurement, reference-frame model
/// vector, and the 1-sigma angular measurement noise (rad).
#[derive(Clone, Copy, Debug)]
pub struct Obs {
    pub b: V3,
    pub r: V3,
    pub sigma: f64,
}

/// Wahba weights a_i proportional to 1/sigma_i^2, normalised so sum a_i = 1.
/// With that normalisation lambda_max is close to 1 and the QUEST Newton
/// iteration has a well-scaled starting point.
pub fn weights(obs: &[Obs]) -> Vec<f64> {
    let inv: Vec<f64> = obs.iter().map(|o| 1.0 / (o.sigma * o.sigma)).collect();
    let s: f64 = inv.iter().sum();
    inv.iter().map(|w| w / s).collect()
}

pub fn attitude_profile(obs: &[Obs], a: &[f64]) -> M3 {
    let mut b = [[0.0; 3]; 3];
    for (o, &w) in obs.iter().zip(a) {
        b = mat_add(&b, &mat_scale(&outer(o.b, o.r), w));
    }
    b
}

/// Wahba loss for a candidate attitude.
pub fn wahba_loss(obs: &[Obs], a_mat: &M3) -> f64 {
    let a = weights(obs);
    obs.iter()
        .zip(&a)
        .map(|(o, &w)| {
            let d = add(o.b, scale(mat_vec(a_mat, o.r), -1.0));
            0.5 * w * dot(d, d)
        })
        .sum()
}

/// TRIAD (Black 1964). Uses only two observations; the first one is trusted
/// fully and the second only fixes the rotation about it, so the more
/// accurate sensor should go first.
pub fn triad(b1: V3, b2: V3, r1: V3, r2: V3) -> M3 {
    let frame = |v1: V3, v2: V3| {
        let t1 = normalize(v1);
        let t2 = normalize(cross(v1, v2));
        let t3 = cross(t1, t2);
        [t1, t2, t3]
    };
    let tb = frame(b1, b2);
    let tr = frame(r1, r2);
    let mut a = [[0.0; 3]; 3];
    for k in 0..3 {
        a = mat_add(&a, &outer(tb[k], tr[k]));
    }
    a
}

/// Davenport's K matrix. The optimal quaternion is the eigenvector of K
/// belonging to its largest eigenvalue.
pub fn davenport_k(b: &M3) -> [[f64; 4]; 4] {
    let s = mat_add(b, &transpose(b));
    let sigma = trace(b);
    let z = [b[1][2] - b[2][1], b[2][0] - b[0][2], b[0][1] - b[1][0]];
    let mut k = [[0.0; 4]; 4];
    for i in 0..3 {
        for j in 0..3 {
            k[i][j] = s[i][j] - if i == j { sigma } else { 0.0 };
        }
        k[i][3] = z[i];
        k[3][i] = z[i];
    }
    k[3][3] = sigma;
    k
}

/// Cyclic Jacobi eigen-decomposition of a symmetric 4x4 matrix.
/// Returns (eigenvalues, eigenvectors as columns).
pub fn jacobi_eig4(m: [[f64; 4]; 4]) -> ([f64; 4], [[f64; 4]; 4]) {
    let mut a = m;
    let mut v = [[0.0; 4]; 4];
    for i in 0..4 {
        v[i][i] = 1.0;
    }
    for _sweep in 0..50 {
        let off: f64 = (0..4)
            .flat_map(|i| (0..4).filter(move |&j| j != i).map(move |j| (i, j)))
            .map(|(i, j)| a[i][j] * a[i][j])
            .sum();
        if off < 1e-30 {
            break;
        }
        for p in 0..3 {
            for q in (p + 1)..4 {
                if a[p][q].abs() < 1e-300 {
                    continue;
                }
                // Rotation angle that zeros a[p][q] (Numerical Recipes form).
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..4 {
                    let akp = a[k][p];
                    let akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..4 {
                    let apk = a[p][k];
                    let aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for k in 0..4 {
                    let vkp = v[k][p];
                    let vkq = v[k][q];
                    v[k][p] = c * vkp - s * vkq;
                    v[k][q] = s * vkp + c * vkq;
                }
            }
        }
    }
    ([a[0][0], a[1][1], a[2][2], a[3][3]], v)
}

/// Davenport q-method: full symmetric eigen-solve of K. The reference answer.
pub fn q_method(obs: &[Obs]) -> (Quat, f64) {
    let a = weights(obs);
    let k = davenport_k(&attitude_profile(obs, &a));
    let (vals, vecs) = jacobi_eig4(k);
    let imax = (0..4).max_by(|&i, &j| vals[i].partial_cmp(&vals[j]).unwrap()).unwrap();
    let q = [vecs[0][imax], vecs[1][imax], vecs[2][imax], vecs[3][imax]];
    (quat_normalize(q), vals[imax])
}

/// Result of a QUEST solve, with diagnostics.
#[derive(Clone, Copy, Debug)]
pub struct QuestResult {
    pub q: Quat,
    pub lambda_max: f64,
    pub newton_iters: usize,
    /// Axis of the 180 deg reference-frame rotation used, if any.
    pub rotated_axis: Option<usize>,
}

/// One QUEST pass on a given profile matrix (no singularity handling).
/// Returns the un-normalised [X, gamma] and the Newton data.
fn quest_core(b: &M3, lambda0: f64) -> ([f64; 4], f64, usize) {
    let s = mat_add(b, &transpose(b));
    let sigma = trace(b);
    let z = [b[1][2] - b[2][1], b[2][0] - b[0][2], b[0][1] - b[1][0]];
    let kappa = trace(&adj3(&s));
    let delta = det3(&s);
    let sz = mat_vec(&s, z);
    let ssz = mat_vec(&s, sz);

    // Characteristic polynomial of K (Shuster & Oh 1981):
    //   psi(l) = l^4 - (a + b) l^2 - c l + (a b + c sigma - d) = 0
    let ca = sigma * sigma - kappa;
    let cb = sigma * sigma + dot(z, z);
    let cc = delta + dot(z, sz);
    let cd = dot(z, ssz);
    let c0 = ca * cb + cc * sigma - cd;

    // Newton from lambda0 = sum a_i; converges in 1-3 steps because
    // lambda_max differs from lambda0 only by the (small) Wahba loss.
    let mut l = lambda0;
    let mut iters = 0;
    for _ in 0..20 {
        let f = l.powi(4) - (ca + cb) * l * l - cc * l + c0;
        let fp = 4.0 * l.powi(3) - 2.0 * (ca + cb) * l - cc;
        let step = f / fp;
        l -= step;
        iters += 1;
        if step.abs() < 1e-14 * lambda0 {
            break;
        }
    }

    // Optimal Gibbs vector via the adjugate of (omega I - S), written so no
    // division happens before the final normalisation.
    //   X = (alpha I + beta S + S^2) z,  gamma = det(omega I - S)
    // with omega = lambda + sigma; by Cayley-Hamilton, using tr S = 2 sigma:
    let alpha = l * l - sigma * sigma + kappa;
    let beta = l - sigma;
    let gamma = (l + sigma) * alpha - delta;
    let x = add(add(scale(z, alpha), scale(sz, beta)), ssz);
    ([x[0], x[1], x[2], gamma], l, iters)
}

/// QUEST (Shuster 1978). Solves for lambda_max with Newton on the
/// characteristic polynomial instead of a full eigen-decomposition.
///
/// The Gibbs-vector form is singular for a 180 deg rotation (gamma -> 0,
/// and at exactly 180 deg X -> 0 too, so the normalised result is pure
/// round-off). Sequential rotations fix this: solve the problem with the
/// reference vectors rotated 180 deg about a coordinate axis, then undo it.
///
/// Which frame to solve in is chosen before solving. Rotating the
/// references 180 deg about axis i turns tr(B) into 2 B_ii - tr(B), which
/// is exactly K_ii, so the four candidate frames correspond to the four
/// diagonal entries of K. The largest one is the frame in which the
/// attitude is closest to identity (largest |q4|). The normalised gamma is
/// still checked, and the next frame tried, for strongly anisotropic
/// observation sets where that heuristic can mislead.
pub fn quest(obs: &[Obs]) -> QuestResult {
    let a = weights(obs);
    let lambda0: f64 = a.iter().sum();
    let b = attitude_profile(obs, &a);
    let tr = trace(&b);
    let kdiag = [2.0 * b[0][0] - tr, 2.0 * b[1][1] - tr, 2.0 * b[2][2] - tr, tr];
    let mut order = [3usize, 0, 1, 2];
    order.sort_by(|&i, &j| kdiag[j].partial_cmp(&kdiag[i]).unwrap());

    let mut best: Option<(f64, QuestResult)> = None;
    for &cand in &order {
        // 180 deg about axis i: r' = R_i r with R_i = diag(+-1), so
        // B' = B R_i^T = B R_i.
        let mut ri = identity();
        if cand < 3 {
            for j in 0..3 {
                ri[j][j] = if j == cand { 1.0 } else { -1.0 };
            }
        }
        let (h, l, it) = quest_core(&mat_mul(&b, &ri), lambda0);
        let hn = (h[0] * h[0] + h[1] * h[1] + h[2] * h[2] + h[3] * h[3]).sqrt();
        let g = (h[3] / hn).abs();
        // b = A' r' = A' R_i r  =>  A = A' R_i
        let q = if cand < 3 {
            dcm_to_quat(&mat_mul(&quat_to_dcm(quat_normalize(h)), &ri))
        } else {
            quat_normalize(h)
        };
        let res = QuestResult {
            q,
            lambda_max: l,
            newton_iters: it,
            rotated_axis: if cand < 3 { Some(cand) } else { None },
        };
        if g > 0.1 {
            return res;
        }
        if best.map_or(true, |(bg, _)| g > bg) {
            best = Some((g, res));
        }
    }
    best.unwrap().1
}

/// QUEST measurement model covariance of the body-frame error vector:
///   P = [ sum_i sigma_i^-2 (I - b_i b_i^T) ]^-1
/// This is the Cramer-Rao bound for the noise model used here, so the
/// q-method and QUEST should both reach it.
pub fn quest_covariance(obs: &[Obs]) -> M3 {
    let mut f = [[0.0; 3]; 3];
    for o in obs {
        let b = normalize(o.b);
        let proj = mat_add(&identity(), &mat_scale(&outer(b, b), -1.0));
        f = mat_add(&f, &mat_scale(&proj, 1.0 / (o.sigma * o.sigma)));
    }
    inv3(&f)
}

/// TRIAD covariance (Shuster & Oh 1981, eq. 37) with observation 1 as anchor:
///   P = s1^2 I + |b1 x b2|^-2 [ s1^2 (b1.b2)(b1 b2^T + b2 b1^T) + (s2^2 - s1^2) b1 b1^T ]
pub fn triad_covariance(b1: V3, s1: f64, b2: V3, s2: f64) -> M3 {
    let b1 = normalize(b1);
    let b2 = normalize(b2);
    let c = cross(b1, b2);
    let c2 = dot(c, c);
    let sym = mat_add(&outer(b1, b2), &outer(b2, b1));
    let term = mat_add(
        &mat_scale(&sym, s1 * s1 * dot(b1, b2)),
        &mat_scale(&outer(b1, b1), s2 * s2 - s1 * s1),
    );
    mat_add(&mat_scale(&identity(), s1 * s1), &mat_scale(&term, 1.0 / c2))
}

// ---------- noise generation ----------

/// splitmix64 + Box-Muller. Deterministic so the Monte Carlo is repeatable.
pub struct Rng {
    state: u64,
    spare: Option<f64>,
}

impl Rng {
    pub fn new(seed: u64) -> Self {
        Rng { state: seed, spare: None }
    }

    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in (0, 1).
    pub fn uniform(&mut self) -> f64 {
        ((self.next_u64() >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    }

    pub fn gauss(&mut self) -> f64 {
        if let Some(s) = self.spare.take() {
            return s;
        }
        let u1 = self.uniform();
        let u2 = self.uniform();
        let r = (-2.0 * u1.ln()).sqrt();
        let th = 2.0 * std::f64::consts::PI * u2;
        self.spare = Some(r * th.sin());
        r * th.cos()
    }

    pub fn unit_vector(&mut self) -> V3 {
        normalize([self.gauss(), self.gauss(), self.gauss()])
    }
}

/// Simulated measurement: true body vector plus Gaussian noise of 1-sigma
/// `sigma` (rad) on each of the two axes perpendicular to it.
pub fn measure(rng: &mut Rng, b_true: V3, sigma: f64) -> V3 {
    let n = [rng.gauss(), rng.gauss(), rng.gauss()];
    let n_perp = add(n, scale(b_true, -dot(n, b_true)));
    normalize(add(b_true, scale(n_perp, sigma)))
}
