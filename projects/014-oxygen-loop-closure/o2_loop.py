# Astron Built this AP
# Project 014 | ECLSS | Oxygen loop closure: PEM electrolysis and Sabatier CO2 reduction
# Author: Aditya (aditya-astron)
# Date: 2026-11-04 IST
"""Steady-state O2 loop of a crewed habitat: how much water does the
Sabatier reactor give back to the electrolyser, and what does it cost?

Loop, per day:

    electrolyser   2 H2O(l) -> 2 H2 + O2         n_w = 2 n_O2, n_H2 = 2 n_O2
    crew           O2 -> CO2  (RQ = n_CO2/n_O2 ~ 0.9)
    Sabatier       CO2 + 4 H2 -> CH4 + 2 H2O     extent xi, at chemical equilibrium
    condenser      2 xi H2O returned, CH4 and leftover CO2/H2 vented

Closure = water returned / water electrolysed = 2 xi / (2 n_O2) = xi / n_O2.
The crew's H2/CO2 ratio is 2/RQ ~ 2.2, well under the 4:1 the reaction
wants, so H2 limits: xi <= n_H2 / 4 = n_O2 / 2 and closure can never pass
50%. Half the hydrogen leaves bound in methane. The equilibrium model says
how close a real reactor temperature gets to that ceiling.

Electrolyser cell voltage (PEM, per cell):

    V = E_rev(T) + (R T / (alpha F)) asinh(i / (2 i0)) + i ASR

E_rev from dG of liquid-water splitting. Stack heat is (V - V_tn) I with
V_tn = dH / (2F), the thermoneutral voltage.

Side reactions (reverse water-gas shift to CO) are ignored; below ~400 C
with a Ni or Ru catalyst they take well under 1% of the carbon.
"""
from dataclasses import dataclass
import math

from scipy.optimize import brentq

from thermo import FARADAY, R, ELECTROLYSIS, SABATIER, ln_K, reaction

M_O2, M_CO2, M_H2O, M_CH4 = 31.998e-3, 44.009e-3, 18.015e-3, 16.043e-3   # kg/mol
DAY = 86400.0
P_STD = 1.0e5   # Pa, standard state for K


def sabatier_extent(n_co2, n_h2, T, P=P_STD, n_inert=0.0):
    """Equilibrium extent xi (same units as the feeds) of CO2 + 4H2 = CH4 + 2H2O.

    Kp = y_CH4 y_H2O^2 / (y_CO2 y_H2^4) * (P/P0)^-2. ln Q(xi) is strictly
    increasing on (0, xi_max), so a bracketed root is unique.
    """
    xi_max = min(n_co2, n_h2 / 4.0)
    lnK = ln_K(SABATIER, T)
    lnP = math.log(P / P_STD)

    def f(xi):
        n_tot = n_co2 + n_h2 + n_inert - 2.0 * xi
        lnQ = (math.log(xi) + 2.0 * math.log(2.0 * xi) + 2.0 * math.log(n_tot)
               - math.log(n_co2 - xi) - 4.0 * math.log(n_h2 - 4.0 * xi) - 2.0 * lnP)
        return lnQ - lnK

    eps = 1e-15 * xi_max
    return brentq(f, eps, xi_max * (1.0 - 1e-15), xtol=1e-18 * xi_max, rtol=1e-14, maxiter=500)


@dataclass
class PEMCell:
    T: float = 353.15          # K, stack temperature
    i0: float = 1.0e-7         # A/cm^2, anode exchange current density (Pt/Ir)
    alpha: float = 1.0         # effective anode transfer coeff., Tafel slope ~70 mV/dec at 80 C
    asr: float = 0.18          # ohm cm^2, membrane + contacts (Nafion 117 class)

    def e_rev(self):
        return reaction(ELECTROLYSIS, self.T)[2] / (2.0 * FARADAY)

    def v_tn(self):
        return reaction(ELECTROLYSIS, self.T)[0] / (2.0 * FARADAY)

    def voltage(self, i):
        act = R * self.T / (self.alpha * FARADAY) * math.asinh(i / (2.0 * self.i0))
        return self.e_rev() + act + i * self.asr


@dataclass
class Crew:
    size: int = 4
    o2_kg_day: float = 0.84     # per crew member, NASA BVAD nominal
    co2_kg_day: float = 1.04

    def n_o2(self):
        return self.size * self.o2_kg_day / M_O2      # mol/day

    def n_co2(self):
        return self.size * self.co2_kg_day / M_CO2


def loop_balance(crew, T_sab, P_sab=P_STD, cell=None, n_cells=12, area_cm2=40.0):
    """Daily mass and energy balance of the closed O2 loop."""
    cell = cell or PEMCell()
    n_o2 = crew.n_o2()
    n_h2 = 2.0 * n_o2
    n_w = 2.0 * n_o2
    n_co2 = crew.n_co2()

    xi = sabatier_extent(n_co2, n_h2, T_sab, P_sab)
    w_back = 2.0 * xi

    # Faraday: 4 electrons per O2. Cells are in series, so each carries the
    # stack current and contributes its share of the total charge.
    I_total = 4.0 * FARADAY * n_o2 / DAY            # A summed over cells
    I_stack = I_total / n_cells
    i_cell = I_stack / area_cm2
    v = cell.voltage(i_cell)
    p_stack = v * I_total
    q_stack = (v - cell.v_tn()) * I_total
    q_sab = -reaction(SABATIER, T_sab)[0] * xi / DAY     # W released

    return {
        "o2_kg": n_o2 * M_O2,
        "water_electrolysed_kg": n_w * M_H2O,
        "water_recovered_kg": w_back * M_H2O,
        "makeup_water_kg": (n_w - w_back) * M_H2O,
        "closure": w_back / n_w,
        "h2_conversion": 4.0 * xi / n_h2,
        "co2_reduced_kg": xi * M_CO2,
        "co2_vented_kg": (n_co2 - xi) * M_CO2,
        "ch4_vented_kg": xi * M_CH4,
        "h2_co2_ratio": n_h2 / n_co2,
        "stack_current_A": I_stack,
        "current_density": i_cell,
        "cell_voltage": v,
        "stack_power_W": p_stack,
        "stack_heat_W": q_stack,
        "sabatier_heat_W": q_sab,
        "kWh_per_kg_o2": p_stack * DAY / 3.6e6 / (n_o2 * M_O2),
    }


def main():
    crew = Crew()
    cell = PEMCell()
    print(f"Crew {crew.size}: O2 {crew.size * crew.o2_kg_day:.2f} kg/d, "
          f"CO2 {crew.size * crew.co2_kg_day:.2f} kg/d, RQ {crew.n_co2() / crew.n_o2():.3f}")
    print(f"PEM cell at {cell.T - 273.15:.0f} C: E_rev {cell.e_rev():.4f} V, "
          f"V_tn {cell.v_tn():.4f} V\n")

    print("Cell polarisation")
    print("  i A/cm2   V cell   eff(HHV)")
    for i in (0.1, 0.5, 1.0, 1.5, 2.0):
        v = cell.voltage(i)
        print(f"  {i:7.2f}  {v:7.4f}  {cell.v_tn() / v:8.3f}")

    print("\nStoichiometric 4:1 feed, 1 bar: equilibrium CO2 conversion")
    for Tc in (200, 300, 400, 500, 600, 700):
        T = Tc + 273.15
        xi = sabatier_extent(1.0, 4.0, T)
        print(f"  {Tc:4d} C  ln K {ln_K(SABATIER, T):7.2f}  X_CO2 {xi:.4f}")

    print("\nClosed loop, crew-limited H2 feed, 1 bar")
    print("  T_sab C  H2 conv  closure  makeup kg/d  CH4 kg/d  Sab heat W")
    for Tc in (250, 300, 350, 400, 450, 500, 600):
        r = loop_balance(crew, Tc + 273.15, cell=cell)
        print(f"  {Tc:7d}  {r['h2_conversion']:7.4f}  {r['closure']:7.4f}  "
              f"{r['makeup_water_kg']:11.3f}  {r['ch4_vented_kg']:8.3f}  {r['sabatier_heat_W']:10.1f}")

    r = loop_balance(crew, 573.15, cell=cell)
    print(f"\nDesign point: Sabatier 300 C, PEM stack 12 cells x 40 cm2 at 80 C")
    print(f"  H2/CO2 feed ratio      {r['h2_co2_ratio']:.3f}  (4.0 needed for full CO2 reduction)")
    print(f"  water electrolysed     {r['water_electrolysed_kg']:.3f} kg/d")
    print(f"  water recovered        {r['water_recovered_kg']:.3f} kg/d  ({100 * r['closure']:.2f}%, ceiling 50%)")
    print(f"  makeup water           {r['makeup_water_kg']:.3f} kg/d  "
          f"= {r['makeup_water_kg'] * 365:.0f} kg/yr")
    print(f"  CO2 vented             {r['co2_vented_kg']:.3f} kg/d")
    print(f"  stack current          {r['stack_current_A']:.2f} A  ({r['current_density']:.3f} A/cm2, "
          f"{r['cell_voltage']:.3f} V/cell)")
    print(f"  stack power            {r['stack_power_W']:.1f} W  ({r['kWh_per_kg_o2']:.2f} kWh/kg O2)")
    print(f"  stack waste heat       {r['stack_heat_W']:.1f} W")
    print(f"  Sabatier heat release  {r['sabatier_heat_W']:.1f} W")


if __name__ == "__main__":
    main()
