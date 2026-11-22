// Astron Built this AP
// Project 056 | Analog habitat systems / ECLSS | Multi-zone CO2 ventilation network for a sealed habitat
// Author: Aditya (aditya-astron)
// Date: 2026-11-22 IST
//
// Small dense linear algebra, enough for a handful of zones: LU solve,
// matrix exponential and a Jacobi symmetric eigensolver. Kept dependency-free
// so the study builds offline.

#[derive(Clone, Debug, PartialEq)]
pub struct Mat {
    pub n: usize,
    pub a: Vec<f64>, // row-major n x n
}

impl Mat {
    pub fn zeros(n: usize) -> Self {
        Mat { n, a: vec![0.0; n * n] }
    }

    pub fn eye(n: usize) -> Self {
        let mut m = Mat::zeros(n);
        for i in 0..n {
            m[(i, i)] = 1.0;
        }
        m
    }

    pub fn mul(&self, b: &Mat) -> Mat {
        let n = self.n;
        let mut c = Mat::zeros(n);
        for i in 0..n {
            for k in 0..n {
                let aik = self[(i, k)];
                if aik == 0.0 {
                    continue;
                }
                for j in 0..n {
                    c.a[i * n + j] += aik * b.a[k * n + j];
                }
            }
        }
        c
    }

    pub fn mul_vec(&self, x: &[f64]) -> Vec<f64> {
        (0..self.n)
            .map(|i| (0..self.n).map(|j| self[(i, j)] * x[j]).sum())
            .collect()
    }

    pub fn scale(&self, s: f64) -> Mat {
        Mat { n: self.n, a: self.a.iter().map(|v| v * s).collect() }
    }

    pub fn add(&self, b: &Mat) -> Mat {
        Mat { n: self.n, a: self.a.iter().zip(&b.a).map(|(x, y)| x + y).collect() }
    }

    pub fn norm1(&self) -> f64 {
        (0..self.n)
            .map(|j| (0..self.n).map(|i| self[(i, j)].abs()).sum::<f64>())
            .fold(0.0, f64::max)
    }

    pub fn trace(&self) -> f64 {
        (0..self.n).map(|i| self[(i, i)]).sum()
    }
}

impl std::ops::Index<(usize, usize)> for Mat {
    type Output = f64;
    fn index(&self, (i, j): (usize, usize)) -> &f64 {
        &self.a[i * self.n + j]
    }
}

impl std::ops::IndexMut<(usize, usize)> for Mat {
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut f64 {
        &mut self.a[i * self.n + j]
    }
}

/// Solve A x = b by Gaussian elimination with partial pivoting.
pub fn solve(a: &Mat, b: &[f64]) -> Vec<f64> {
    let n = a.n;
    let mut m = a.a.clone();
    let mut x = b.to_vec();
    for k in 0..n {
        let p = (k..n)
            .max_by(|&i, &j| m[i * n + k].abs().partial_cmp(&m[j * n + k].abs()).unwrap())
            .unwrap();
        assert!(m[p * n + k].abs() > 1e-300, "singular matrix");
        if p != k {
            for j in 0..n {
                m.swap(k * n + j, p * n + j);
            }
            x.swap(k, p);
        }
        for i in k + 1..n {
            let f = m[i * n + k] / m[k * n + k];
            if f == 0.0 {
                continue;
            }
            for j in k..n {
                m[i * n + j] -= f * m[k * n + j];
            }
            x[i] -= f * x[k];
        }
    }
    for k in (0..n).rev() {
        let s: f64 = (k + 1..n).map(|j| m[k * n + j] * x[j]).sum();
        x[k] = (x[k] - s) / m[k * n + k];
    }
    x
}

/// exp(A) by scaling and squaring with a truncated Taylor series.
/// After scaling ||A/2^s||_1 <= 1/2, so 20 terms leave a truncation error
/// below 0.5^21/21! ~ 1e-26, far under the squaring round-off.
pub fn expm(a: &Mat) -> Mat {
    let nrm = a.norm1();
    let s = if nrm > 0.5 { (nrm / 0.5).log2().ceil() as i32 } else { 0 };
    let x = a.scale(0.5f64.powi(s));
    let mut term = Mat::eye(a.n);
    let mut sum = Mat::eye(a.n);
    for k in 1..=20 {
        term = term.mul(&x).scale(1.0 / k as f64);
        sum = sum.add(&term);
    }
    for _ in 0..s {
        sum = sum.mul(&sum);
    }
    sum
}

/// Eigenvalues of a real symmetric matrix by cyclic Jacobi rotations,
/// returned in ascending order.
pub fn sym_eigenvalues(s: &Mat) -> Vec<f64> {
    let n = s.n;
    let mut a = s.clone();
    for _sweep in 0..100 {
        let off: f64 = (0..n)
            .flat_map(|i| (0..n).filter(move |&j| j != i).map(move |j| (i, j)))
            .map(|(i, j)| a[(i, j)] * a[(i, j)])
            .sum();
        if off < 1e-30 * a.norm1().powi(2) {
            break;
        }
        for p in 0..n {
            for q in p + 1..n {
                if a[(p, q)].abs() < 1e-300 {
                    continue;
                }
                // Rotation angle that zeroes a_pq (Golub & Van Loan 8.5.2).
                let tau = (a[(q, q)] - a[(p, p)]) / (2.0 * a[(p, q)]);
                let t = tau.signum() / (tau.abs() + (1.0 + tau * tau).sqrt());
                let t = if tau == 0.0 { 1.0 } else { t };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let sn = t * c;
                for k in 0..n {
                    let akp = a[(k, p)];
                    let akq = a[(k, q)];
                    a[(k, p)] = c * akp - sn * akq;
                    a[(k, q)] = sn * akp + c * akq;
                }
                for k in 0..n {
                    let apk = a[(p, k)];
                    let aqk = a[(q, k)];
                    a[(p, k)] = c * apk - sn * aqk;
                    a[(q, k)] = sn * apk + c * aqk;
                }
            }
        }
    }
    let mut ev: Vec<f64> = (0..n).map(|i| a[(i, i)]).collect();
    ev.sort_by(|x, y| x.partial_cmp(y).unwrap());
    ev
}
