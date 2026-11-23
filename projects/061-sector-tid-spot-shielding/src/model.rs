// Astron Built this AP
// Project 061 | Radiation shielding | Sector-based TID analysis and spot-shield sizing for an avionics box
// Author: Aditya (aditya-astron)
// Date: 2026-11-23 IST
//
// Mass model of a 1.6 m GEO bus with an avionics box bolted to the inside of
// the +X panel. Origin at the bus centre, lengths in cm, densities in g/cm^3
// (aluminium-equivalent by mass thickness). Component densities are smeared
// averages: propellant at a mission-average 50 % fill, battery cells and
// voids, a reaction-wheel and harness block.

use crate::dose::AL_RHO;
use crate::geometry::{Scene, V3};

pub const PANEL_HALF: f64 = 80.0;
pub const PANEL_WALL: f64 = 1.0;
pub const PANEL_RHO: f64 = 0.40; // honeycomb + facesheets: 0.40 g/cm^2

pub const BOX_LO: V3 = [49.0, -10.0, -7.5];
pub const BOX_HI: V3 = [79.0, 10.0, 7.5];
pub const BOX_WALL: f64 = 0.20; // 2 mm Al

pub const CARD_X: [f64; 5] = [52.0, 57.0, 62.0, 67.0, 72.0];
pub const CARD_T: f64 = 0.16;
pub const CARD_RHO: f64 = 1.85; // FR-4 with copper planes
pub const CARD_HALF_Y: f64 = 9.0;
pub const CARD_HALF_Z: f64 = 6.5;

/// Spot-shield can: inner footprint half-width and inner height.
pub const CAN_HALF: f64 = 1.2;
pub const CAN_HEIGHT: f64 = 0.8;

#[derive(Clone, Debug)]
pub struct Part {
    pub name: &'static str,
    pub card: usize,
    pub y: f64,
    pub z: f64,
    /// Part TID capability, krad(Si).
    pub capability: f64,
}

impl Part {
    /// Sensitive volume 1 mm above the card's component side (+X face).
    pub fn point(&self) -> V3 {
        [CARD_X[self.card] + CARD_T + 0.1, self.y, self.z]
    }
}

pub fn parts() -> Vec<Part> {
    vec![
        Part { name: "U1 FPGA (card 5, corner)", card: 4, y: 6.0, z: 4.0, capability: 100.0 },
        Part { name: "U2 SRAM (card 5, centre)", card: 4, y: 0.0, z: 0.0, capability: 10.0 },
        Part { name: "U3 ADC  (card 3, centre)", card: 2, y: 0.0, z: 0.0, capability: 20.0 },
        Part { name: "U4 LDO  (card 1, corner)", card: 0, y: -6.0, z: -3.0, capability: 10.0 },
    ]
}

/// Everything outside the avionics box.
pub fn bus_scene() -> Scene {
    let mut s = Scene::new();
    let h = PANEL_HALF;
    s.add_box_shell([-h, -h, -h], [h, h, h], PANEL_WALL, PANEL_RHO);
    // Propellant tank: 1.5 mm titanium shell, hydrazine at 50 % average fill.
    s.add_sphere_shell([0.0, 0.0, 0.0], 35.0, 0.15, 4.43);
    s.add_sphere([0.0, 0.0, 0.0], 34.85, 0.50);
    // Battery pack and a wheel/harness block either side of the box.
    s.add_box([20.0, 30.0, -20.0], [45.0, 60.0, 20.0], 1.20);
    s.add_box([40.0, -60.0, -15.0], [70.0, -30.0, 15.0], 0.60);
    s
}

/// Avionics box with all six walls thickened inward by `extra` (cm Al) and
/// an aluminium doubler of thickness `doubler` on the inside of the +X wall.
pub fn avionics_box(extra: f64, doubler: f64) -> Scene {
    let mut s = Scene::new();
    s.add_box_shell(BOX_LO, BOX_HI, BOX_WALL + extra, AL_RHO);
    if doubler > 0.0 {
        let x1 = BOX_HI[0] - BOX_WALL - extra;
        let w = BOX_WALL + extra;
        s.add_box([x1 - doubler, BOX_LO[1] + w, BOX_LO[2] + w], [x1, BOX_HI[1] - w, BOX_HI[2] - w], AL_RHO);
    }
    for &x in CARD_X.iter() {
        s.add_box([x, -CARD_HALF_Y, -CARD_HALF_Z], [x + CARD_T, CARD_HALF_Y, CARD_HALF_Z], CARD_RHO);
    }
    s
}

/// Aluminium spot-shield can of wall w over a part: lid plus four side
/// walls standing on the card. With `backside`, a plate of the same
/// thickness and footprint is added under the card.
pub fn spot_can(part: &Part, w: f64, backside: bool) -> Scene {
    let mut s = Scene::new();
    if w <= 0.0 {
        return s;
    }
    let top = CARD_X[part.card] + CARD_T;
    let (y, z, c, h) = (part.y, part.z, CAN_HALF, CAN_HEIGHT);
    s.add_box([top + h, y - c - w, z - c - w], [top + h + w, y + c + w, z + c + w], AL_RHO);
    s.add_box([top, y + c, z - c - w], [top + h, y + c + w, z + c + w], AL_RHO);
    s.add_box([top, y - c - w, z - c - w], [top + h, y - c, z + c + w], AL_RHO);
    s.add_box([top, y - c, z + c], [top + h, y + c, z + c + w], AL_RHO);
    s.add_box([top, y - c, z - c - w], [top + h, y + c, z - c], AL_RHO);
    if backside {
        let bot = CARD_X[part.card];
        s.add_box([bot - w, y - c - w, z - c - w], [bot, y + c + w, z + c + w], AL_RHO);
    }
    s
}

/// Index of the avionics-box face a ray from p (inside the box) leaves
/// through: 0..5 = +X, -X, +Y, -Y, +Z, -Z.
pub fn exit_face(p: &V3, u: &V3) -> usize {
    let mut best = (f64::INFINITY, 0);
    for k in 0..3 {
        if u[k] == 0.0 {
            continue;
        }
        let (bound, face) = if u[k] > 0.0 { (BOX_HI[k], 2 * k) } else { (BOX_LO[k], 2 * k + 1) };
        let s = (bound - p[k]) / u[k];
        if s < best.0 {
            best = (s, face);
        }
    }
    best.1
}
