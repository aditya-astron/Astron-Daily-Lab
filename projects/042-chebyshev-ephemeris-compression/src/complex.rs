// Astron Built this AP
// Project 042 | Numerical methods / Navigation | Chebyshev ephemeris compression for an eccentric orbit
// Author: Aditya (aditya-astron)
// Date: 2026-11-16 IST

//! Minimal complex type, enough for locating Bernstein ellipses.

use std::ops::{Add, Mul, Sub};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct C64 {
    pub re: f64,
    pub im: f64,
}

impl C64 {
    pub fn new(re: f64, im: f64) -> Self {
        C64 { re, im }
    }
    pub fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }
    /// Principal square root, computed without cancellation.
    pub fn sqrt(self) -> Self {
        let r = self.abs();
        if r == 0.0 {
            return self;
        }
        let t = (0.5 * (r + self.re.abs())).sqrt();
        if self.re >= 0.0 {
            C64::new(t, self.im / (2.0 * t))
        } else {
            C64::new(self.im.abs() / (2.0 * t), t.copysign(self.im))
        }
    }
}

impl Add for C64 {
    type Output = C64;
    fn add(self, o: C64) -> C64 {
        C64::new(self.re + o.re, self.im + o.im)
    }
}
impl Sub for C64 {
    type Output = C64;
    fn sub(self, o: C64) -> C64 {
        C64::new(self.re - o.re, self.im - o.im)
    }
}
impl Mul for C64 {
    type Output = C64;
    fn mul(self, o: C64) -> C64 {
        C64::new(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
    }
}
