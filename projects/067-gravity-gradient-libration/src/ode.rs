// Astron Built this AP
// Project 067 | Attitude control | Gravity-gradient stability regions and eccentric-orbit pitch libration
// Author: Aditya (aditya-astron)
// Date: 2026-11-27 IST
//
// Fixed-step classical RK4. Every problem here is smooth and periodic, so a
// fixed step lets the monodromy integrals land exactly on 2*pi and makes the
// step-halving convergence checks in the tests clean fourth order.

pub fn rk4_step<F>(f: &F, t: f64, y: &mut [f64], h: f64)
where
    F: Fn(f64, &[f64], &mut [f64]),
{
    let n = y.len();
    let mut k1 = vec![0.0; n];
    let mut k2 = vec![0.0; n];
    let mut k3 = vec![0.0; n];
    let mut k4 = vec![0.0; n];
    let mut tmp = vec![0.0; n];

    f(t, y, &mut k1);
    for i in 0..n {
        tmp[i] = y[i] + 0.5 * h * k1[i];
    }
    f(t + 0.5 * h, &tmp, &mut k2);
    for i in 0..n {
        tmp[i] = y[i] + 0.5 * h * k2[i];
    }
    f(t + 0.5 * h, &tmp, &mut k3);
    for i in 0..n {
        tmp[i] = y[i] + h * k3[i];
    }
    f(t + h, &tmp, &mut k4);
    for i in 0..n {
        y[i] += h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
}

/// Integrate from t0 to t1 in `steps` equal steps, calling `obs` after each.
pub fn rk4<F, O>(f: &F, t0: f64, t1: f64, y: &mut [f64], steps: usize, mut obs: O)
where
    F: Fn(f64, &[f64], &mut [f64]),
    O: FnMut(f64, &[f64]),
{
    let h = (t1 - t0) / steps as f64;
    for i in 0..steps {
        let t = t0 + i as f64 * h;
        rk4_step(f, t, y, h);
        obs(t + h, y);
    }
}

/// State transition matrix of a linear system x' = A(t) x over [t0, t1],
/// column by column. Returned row-major, n x n.
pub fn transition_matrix<F>(f: &F, n: usize, t0: f64, t1: f64, steps: usize) -> Vec<f64>
where
    F: Fn(f64, &[f64], &mut [f64]),
{
    let mut phi = vec![0.0; n * n];
    for j in 0..n {
        let mut x = vec![0.0; n];
        x[j] = 1.0;
        rk4(f, t0, t1, &mut x, steps, |_, _| {});
        for i in 0..n {
            phi[i * n + j] = x[i];
        }
    }
    phi
}

/// Characteristic polynomial coefficients c[0..=n] of a square matrix
/// (det(lambda I - M) = sum c[k] lambda^(n-k), c[0] = 1) by Faddeev-LeVerrier.
pub fn char_poly(m: &[f64], n: usize) -> Vec<f64> {
    let mut c = vec![0.0; n + 1];
    c[0] = 1.0;
    let mut mk = vec![0.0; n * n]; // M_k, starts as zero so M_1 = I
    for k in 1..=n {
        // M_k = M * M_{k-1} + c_{k-1} I
        let mut next = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                let mut s = 0.0;
                for l in 0..n {
                    s += m[i * n + l] * mk[l * n + j];
                }
                next[i * n + j] = s;
            }
            next[i * n + i] += c[k - 1];
        }
        mk = next;
        // c_k = -tr(M M_k) / k
        let mut tr = 0.0;
        for i in 0..n {
            for l in 0..n {
                tr += m[i * n + l] * mk[l * n + i];
            }
        }
        c[k] = -tr / k as f64;
    }
    c
}

/// Multiply out prod (lambda - r_i) for complex roots given as (re, im) pairs.
pub fn poly_from_roots(roots: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let mut p = vec![(1.0, 0.0)];
    for &(rr, ri) in roots {
        let mut q = vec![(0.0, 0.0); p.len() + 1];
        for (k, &(pr, pi)) in p.iter().enumerate() {
            q[k].0 += pr;
            q[k].1 += pi;
            q[k + 1].0 -= pr * rr - pi * ri;
            q[k + 1].1 -= pr * ri + pi * rr;
        }
        p = q;
    }
    p
}
