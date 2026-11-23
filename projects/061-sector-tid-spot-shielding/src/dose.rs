// Astron Built this AP
// Project 061 | Radiation shielding | Sector-based TID analysis and spot-shield sizing for an avionics box
// Author: Aditya (aditya-astron)
// Date: 2026-11-23 IST
//
// Dose-depth curve: total ionising dose at the centre of a solid aluminium
// sphere of radius t (g/cm^2), the form SHIELDOSE-2 and similar codes report.
// It is represented as a sum of exponentials,
//
//     D(t) = sum_k A_k exp(-t / lambda_k)
//
// which fits the steep trapped-electron fall-off plus the slow
// proton/bremsstrahlung tail well, and makes the two-sided infinite-slab
// geometry exact: averaging exp(-w/(lambda |mu|)) over the sphere of
// directions gives the exponential integral E2(w/lambda).

#[derive(Clone, Debug)]
pub struct DoseDepth {
    /// (A_k in krad(Si), lambda_k in g/cm^2)
    pub terms: Vec<(f64, f64)>,
}

/// Aluminium: 1 mm = 0.270 g/cm^2.
pub const AL_RHO: f64 = 2.70;

pub fn mm_al_to_areal(mm: f64) -> f64 {
    mm * 0.1 * AL_RHO
}

pub fn areal_to_mm_al(t: f64) -> f64 {
    t / (0.1 * AL_RHO)
}

impl DoseDepth {
    /// Illustrative three-term curve shaped like a 15-year GEO mission:
    /// ~420 krad behind 1 mm Al, ~2 krad behind 10 mm. Replace with the
    /// program's SHIELDOSE-2 / OMERE output for flight work.
    pub fn geo_15yr_illustrative() -> Self {
        DoseDepth {
            terms: vec![(4000.0, 0.10), (300.0, 0.40), (3.0, 5.0)],
        }
    }

    pub fn eval(&self, t: f64) -> f64 {
        self.terms.iter().map(|(a, l)| a * (-t / l).exp()).sum()
    }

    /// Dose midway between two infinite plates of areal density w each.
    pub fn two_sided_slab(&self, w: f64) -> f64 {
        self.terms.iter().map(|(a, l)| a * e2(w / l)).sum()
    }

    /// Solid-sphere thickness that gives dose d (D is strictly decreasing).
    pub fn equivalent_thickness(&self, d: f64) -> f64 {
        assert!(d > 0.0 && d <= self.eval(0.0));
        let (mut lo, mut hi) = (0.0, 1.0);
        while self.eval(hi) > d {
            hi *= 2.0;
        }
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if self.eval(mid) > d { lo = mid } else { hi = mid }
            if hi - lo < 1e-14 * hi.max(1.0) {
                break;
            }
        }
        0.5 * (lo + hi)
    }
}

const EULER_GAMMA: f64 = 0.577_215_664_901_532_9;

/// Exponential integral E1(x), x > 0. Power series below 1, modified Lentz
/// continued fraction above (Numerical Recipes 6.3 with n = 1).
pub fn e1(x: f64) -> f64 {
    assert!(x > 0.0);
    if x <= 1.0 {
        // E1 = -gamma - ln x - sum_{k>=1} (-x)^k / (k k!)
        let mut sum = 0.0;
        let mut term = 1.0;
        for k in 1..60 {
            term *= -x / k as f64;
            let add = term / k as f64;
            sum += add;
            if add.abs() < 1e-17 * sum.abs().max(1e-300) {
                break;
            }
        }
        -EULER_GAMMA - x.ln() - sum
    } else {
        let tiny = 1e-300;
        let mut b = x + 1.0;
        let mut c = 1.0 / tiny;
        let mut d = 1.0 / b;
        let mut h = d;
        for i in 1..500 {
            let a = -((i * i) as f64);
            b += 2.0;
            d = 1.0 / (a * d + b);
            c = b + a / c;
            let del = c * d;
            h *= del;
            if (del - 1.0).abs() < 1e-16 {
                break;
            }
        }
        h * (-x).exp()
    }
}

/// E2(x) = integral_0^1 exp(-x/mu) dmu = exp(-x) - x E1(x).
pub fn e2(x: f64) -> f64 {
    if x == 0.0 {
        return 1.0;
    }
    (-x).exp() - x * e1(x)
}
