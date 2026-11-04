# Astron Built this AP
# Project 012 | Radiation shielding | Solar proton event dose behind aluminium vs polyethylene shielding
# Author: Aditya (aditya-astron)
# Date: 2026-11-04 IST
"""Proton electronic stopping power (Bethe) and CSDA range tables.

Mass stopping power, MeV cm^2/g:

    S = K (Z/A) / b^2 * [ 0.5 ln(2 me c^2 b^2 g^2 Tmax / I^2) - b^2 ]
    Tmax = 2 me c^2 b^2 g^2 / (1 + 2 g me/M + (me/M)^2)

Shell, Barkas and density corrections are left out. Above ~5 MeV that costs
well under 1% for low-Z materials, which is small next to the spread in
solar-event spectra; below that range contributions are tiny anyway.
"""
from dataclasses import dataclass, field

import numpy as np

K = 0.307075            # 4 pi N_A r_e^2 m_e c^2, MeV cm^2/mol
ME = 0.51099895         # electron rest energy, MeV
MP = 938.27209          # proton rest energy, MeV
E_MIN = 1.0             # bottom of the tables, MeV
E_MAX = 3000.0
BRAGG_KLEEMAN_P = 1.75  # R ~ E^p below E_MIN, so R(E_MIN) = E_MIN / (p S)


def bethe(E, z_over_a, I_eV):
    E = np.asarray(E, dtype=float)
    g = 1.0 + E / MP
    b2 = 1.0 - 1.0 / g**2
    bg2 = b2 * g**2
    r = ME / MP
    tmax = 2 * ME * bg2 / (1 + 2 * g * r + r * r)
    I = I_eV * 1e-6
    L = 0.5 * np.log(2 * ME * bg2 * tmax / I**2) - b2
    return K * z_over_a / b2 * L


@dataclass
class Material:
    """Range-energy tables on a log grid; all lookups are log-log interpolated."""
    name: str
    z_over_a: float
    I_eV: float
    n: int = 4000
    _logE: np.ndarray = field(init=False, repr=False)
    _logR: np.ndarray = field(init=False, repr=False)

    def __post_init__(self):
        E = np.geomspace(E_MIN, E_MAX, self.n)
        S = bethe(E, self.z_over_a, self.I_eV)
        # Range = integral dE/S, trapezoid in ln E (integrand E/S is smooth there).
        f = E / S
        dlnE = np.diff(np.log(E))
        R = np.concatenate(([0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * dlnE)))
        R += E_MIN / (BRAGG_KLEEMAN_P * S[0])
        self._logE, self._logR = np.log(E), np.log(R)

    def S(self, E):
        return bethe(E, self.z_over_a, self.I_eV)

    def range(self, E):
        return np.exp(np.interp(np.log(E), self._logE, self._logR))

    def energy(self, R):
        """Inverse of range(). Ranges above the table return E_MAX (never used in practice)."""
        R = np.maximum(np.asarray(R, dtype=float), 1e-300)
        return np.exp(np.interp(np.log(R), self._logR, self._logE))


# ICRU 37/49 mean excitation energies.
WATER = Material("water", 10 / 18.0153, 75.0)
ALUMINIUM = Material("aluminium", 13 / 26.9815, 166.0)
POLYETHYLENE = Material("polyethylene", 0.57034, 57.4)
