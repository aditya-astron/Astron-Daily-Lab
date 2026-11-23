// Astron Built this AP
// Project 061 | Radiation shielding | Sector-based TID analysis and spot-shield sizing for an avionics box
// Author: Aditya (aditya-astron)
// Date: 2026-11-23 IST
//
// Constructive mass model for sector shielding analysis. Every body is an
// axis-aligned box or a sphere with a density in g/cm^3 (lengths in cm), so
// the chord a ray cuts through it times the density is an areal density in
// g/cm^2. Hollow bodies are an outer solid plus an inner solid of negative
// density: the ray's inner chord cancels exactly the part of the outer chord
// that lies in the cavity. Materials are treated as aluminium-equivalent by
// mass thickness, which is the usual first-cut assumption for TID.

pub type V3 = [f64; 3];

#[derive(Clone, Debug)]
pub enum Solid {
    Box { lo: V3, hi: V3, rho: f64 },
    Sphere { c: V3, r: f64, rho: f64 },
}

/// Length of the half-line p + s*u (s >= 0) inside the box [lo, hi].
/// Slab method: intersect the three parameter intervals.
pub fn chord_box(p: &V3, u: &V3, lo: &V3, hi: &V3) -> f64 {
    let mut s0 = 0.0_f64;
    let mut s1 = f64::INFINITY;
    for k in 0..3 {
        if u[k] == 0.0 {
            if p[k] < lo[k] || p[k] > hi[k] {
                return 0.0;
            }
            continue;
        }
        let inv = 1.0 / u[k];
        let (a, b) = {
            let ta = (lo[k] - p[k]) * inv;
            let tb = (hi[k] - p[k]) * inv;
            if ta < tb { (ta, tb) } else { (tb, ta) }
        };
        s0 = s0.max(a);
        s1 = s1.min(b);
        if s1 <= s0 {
            return 0.0;
        }
    }
    s1 - s0
}

/// Length of the half-line p + s*u (s >= 0, |u| = 1) inside the sphere.
pub fn chord_sphere(p: &V3, u: &V3, c: &V3, r: f64) -> f64 {
    let d = [p[0] - c[0], p[1] - c[1], p[2] - c[2]];
    let b = d[0] * u[0] + d[1] * u[1] + d[2] * u[2];
    let cc = d[0] * d[0] + d[1] * d[1] + d[2] * d[2] - r * r;
    let disc = b * b - cc;
    if disc <= 0.0 {
        return 0.0;
    }
    let q = disc.sqrt();
    let s0 = (-b - q).max(0.0);
    let s1 = -b + q;
    if s1 <= s0 { 0.0 } else { s1 - s0 }
}

#[derive(Clone, Debug, Default)]
pub struct Scene {
    pub solids: Vec<Solid>,
}

impl Scene {
    pub fn new() -> Self {
        Scene { solids: Vec::new() }
    }

    pub fn add_box(&mut self, lo: V3, hi: V3, rho: f64) {
        self.solids.push(Solid::Box { lo, hi, rho });
    }

    /// Closed box shell with outer extent [lo, hi] and uniform wall thickness.
    pub fn add_box_shell(&mut self, lo: V3, hi: V3, wall: f64, rho: f64) {
        self.add_box(lo, hi, rho);
        let ilo = [lo[0] + wall, lo[1] + wall, lo[2] + wall];
        let ihi = [hi[0] - wall, hi[1] - wall, hi[2] - wall];
        self.add_box(ilo, ihi, -rho);
    }

    pub fn add_sphere(&mut self, c: V3, r: f64, rho: f64) {
        self.solids.push(Solid::Sphere { c, r, rho });
    }

    pub fn add_sphere_shell(&mut self, c: V3, r_out: f64, wall: f64, rho: f64) {
        self.add_sphere(c, r_out, rho);
        self.add_sphere(c, r_out - wall, -rho);
    }

    pub fn extend(&mut self, other: &Scene) {
        self.solids.extend(other.solids.iter().cloned());
    }

    /// Areal density (g/cm^2) seen from p along unit direction u.
    pub fn areal_density(&self, p: &V3, u: &V3) -> f64 {
        let mut t = 0.0;
        for s in &self.solids {
            t += match s {
                Solid::Box { lo, hi, rho } => rho * chord_box(p, u, lo, hi),
                Solid::Sphere { c, r, rho } => rho * chord_sphere(p, u, c, *r),
            };
        }
        // Cancellation of an outer and inner chord can leave -1e-15.
        t.max(0.0)
    }

    /// Mass in grams (negative-density cavities subtract).
    pub fn mass(&self) -> f64 {
        self.solids
            .iter()
            .map(|s| match s {
                Solid::Box { lo, hi, rho } => {
                    rho * (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2])
                }
                Solid::Sphere { r, rho, .. } => rho * 4.0 / 3.0 * std::f64::consts::PI * r * r * r,
            })
            .sum()
    }
}
