# Astron Built this AP
# Project 001 | ECLSS | LiOH CO2 Scrubber Sizing for a 4-Crew 30-Day Mission
# Author: Aditya (aditya-astron)
# Date: 2026-10-30 IST
"""
Size a non-regenerable LiOH CO2 removal system for a sealed crew cabin.

Two levels of analysis:
  1. Static stoichiometric budget: how much LiOH the mission consumes.
  2. Dynamic cabin model: ppCO2 vs. time with a single active canister whose
     removal efficiency collapses near breakthrough. A canister is swapped
     when cabin ppCO2 crosses the action limit, so the real usable capacity
     per canister falls out of the simulation instead of being assumed.

Chemistry (overall, anhydrous LiOH):
    2 LiOH(s) + CO2(g) -> Li2CO3(s) + H2O
"""

from dataclasses import dataclass, field
import math

import numpy as np
from scipy.integrate import solve_ivp

R_UNIVERSAL = 8.314462618        # J/(mol K)
MMHG_PER_PA = 1.0 / 133.322368

# Molar masses, kg/mol
M_LIOH = 23.946e-3
M_CO2 = 44.009e-3
M_H2O = 18.015e-3
M_LI2CO3 = 73.891e-3

# Standard enthalpies of formation at 298 K, kJ/mol (NIST-JANAF)
HF_LIOH = -484.9
HF_CO2 = -393.5
HF_LI2CO3 = -1216.0
HF_H2O_GAS = -241.8
HF_H2O_LIQ = -285.8


def lioh_per_kg_co2():
    """Theoretical kg LiOH consumed per kg CO2 absorbed (2:1 molar)."""
    return 2.0 * M_LIOH / M_CO2


def water_per_kg_co2():
    return M_H2O / M_CO2


def reaction_heat_per_kg_co2(water_liquid=False):
    """Heat released (J) per kg CO2 absorbed. Positive = exothermic.

    Liquid product water is the upper bound; in a flow-through bed most of
    it leaves as vapour, which is the number the cabin heat load should use.
    """
    hf_water = HF_H2O_LIQ if water_liquid else HF_H2O_GAS
    dh_rxn = (HF_LI2CO3 + hf_water) - (2.0 * HF_LIOH + HF_CO2)   # kJ/mol CO2
    return -dh_rxn * 1e3 / M_CO2


@dataclass
class Mission:
    crew: int = 4
    days: float = 30.0
    co2_per_crew_day: float = 1.04     # kg, NASA BVAD nominal metabolic rate
    cabin_volume: float = 20.0         # m^3 free gas volume
    cabin_temp: float = 295.15         # K
    ppco2_start_mmhg: float = 0.5
    ppco2_limit_mmhg: float = 3.0      # swap trigger, long-duration exposure limit

    @property
    def co2_rate(self):
        """Crew CO2 generation, kg/s."""
        return self.crew * self.co2_per_crew_day / 86400.0


@dataclass
class Canister:
    lioh_mass: float = 4.0         # kg LiOH per canister
    flow: float = 0.010            # m^3/s cabin air through the bed
    eta_fresh: float = 0.98        # single-pass removal efficiency, fresh bed
    u50: float = 0.88              # utilisation at which efficiency halves
    width: float = 0.03            # breakthrough sharpness (utilisation units)

    @property
    def co2_capacity(self):
        """Theoretical CO2 a canister can hold, kg."""
        return self.lioh_mass / lioh_per_kg_co2()

    def efficiency(self, u):
        """Logistic breakthrough curve, normalised so eta(0) = eta_fresh."""
        f0 = 1.0 / (1.0 + math.exp(-self.u50 / self.width))
        return self.eta_fresh / (1.0 + np.exp((u - self.u50) / self.width)) / f0


def pa_to_kg_m3(p_pa, temp):
    return p_pa * M_CO2 / (R_UNIVERSAL * temp)


def kg_m3_to_mmhg(c, temp):
    return c * R_UNIVERSAL * temp / M_CO2 * MMHG_PER_PA


def static_budget(m: Mission, utilisation=1.0):
    """Stoichiometric LiOH budget assuming all generated CO2 is absorbed."""
    co2_total = m.co2_rate * m.days * 86400.0
    lioh = co2_total * lioh_per_kg_co2() / utilisation
    return {"co2_kg": co2_total, "lioh_kg": lioh}


def steady_state_ppco2(m: Mission, c: Canister, eta=None):
    """Well-mixed cabin equilibrium: G = eta * Q * C  ->  C = G / (eta Q)."""
    eta = c.eta_fresh if eta is None else eta
    return kg_m3_to_mmhg(m.co2_rate / (eta * c.flow), m.cabin_temp)


@dataclass
class SimResult:
    t: np.ndarray
    ppco2: np.ndarray
    swap_times: list = field(default_factory=list)
    swap_utilisation: list = field(default_factory=list)
    canisters_used: int = 0
    co2_generated: float = 0.0
    co2_absorbed: float = 0.0
    co2_cabin_delta: float = 0.0


def simulate(m: Mission, c: Canister, max_step=600.0):
    """Integrate cabin CO2 mass and canister loading, swapping on the limit.

    State y = [C (kg/m^3 cabin CO2), A (kg CO2 in active canister)]
        dC/dt = (G - eta(u) Q C) / V
        dA/dt = eta(u) Q C,           u = A / capacity
    """
    G, V, Q = m.co2_rate, m.cabin_volume, c.flow
    cap = c.co2_capacity
    c_limit = pa_to_kg_m3(m.ppco2_limit_mmhg / MMHG_PER_PA, m.cabin_temp)
    t_end = m.days * 86400.0

    def rhs(_, y):
        removal = c.efficiency(y[1] / cap) * Q * y[0]
        return [(G - removal) / V, removal]

    def over_limit(_, y):
        return y[0] - c_limit
    over_limit.terminal = True
    over_limit.direction = 1

    c0 = pa_to_kg_m3(m.ppco2_start_mmhg / MMHG_PER_PA, m.cabin_temp)
    t0, y = 0.0, [c0, 0.0]
    res = SimResult(t=np.empty(0), ppco2=np.empty(0))
    ts, cs = [], []
    absorbed = 0.0
    res.canisters_used = 1

    while t0 < t_end:
        sol = solve_ivp(rhs, (t0, t_end), y, method="LSODA", events=over_limit,
                        max_step=max_step, rtol=1e-8, atol=[1e-10, 1e-8])
        if not sol.success:
            raise RuntimeError(sol.message)
        ts.append(sol.t)
        cs.append(sol.y[0])
        t0 = sol.t[-1]
        y_end = sol.y[:, -1]
        absorbed += y_end[1]
        if sol.status == 1:                  # limit hit: swap canister
            res.swap_times.append(t0)
            res.swap_utilisation.append(y_end[1] / cap)
            res.canisters_used += 1
            # Restart exactly on the limit; a fresh bed makes dC/dt < 0, so
            # the upward-crossing event cannot re-fire at the same instant.
            y = [y_end[0], 0.0]
        else:
            break

    res.t = np.concatenate(ts)
    res.ppco2 = kg_m3_to_mmhg(np.concatenate(cs), m.cabin_temp)
    res.co2_generated = G * t_end
    res.co2_absorbed = absorbed
    res.co2_cabin_delta = (cs[-1][-1] - c0) * V
    return res


def report(m: Mission, c: Canister, spares=0.10):
    sim = simulate(m, c)
    static = static_budget(m)
    mean_u = float(np.mean(sim.swap_utilisation)) if sim.swap_utilisation else float("nan")
    n_flight = math.ceil(sim.canisters_used * (1.0 + spares))
    absorbed = sim.co2_absorbed
    heat_w = reaction_heat_per_kg_co2() * absorbed / (m.days * 86400.0)

    lines = [
        f"Mission: {m.crew} crew, {m.days:.0f} days, cabin {m.cabin_volume:.1f} m^3",
        f"CO2 generated                : {static['co2_kg']:8.2f} kg ({m.co2_rate*86400:.2f} kg/day)",
        f"LiOH per kg CO2 (theory)     : {lioh_per_kg_co2():8.4f} kg/kg",
        f"LiOH, 100% utilisation       : {static['lioh_kg']:8.2f} kg",
        f"Fresh-bed steady ppCO2       : {steady_state_ppco2(m, c):8.3f} mmHg",
        f"Swap limit                   : {m.ppco2_limit_mmhg:8.3f} mmHg",
        f"Canister CO2 capacity        : {c.co2_capacity:8.3f} kg ({c.lioh_mass:.1f} kg LiOH)",
        f"Mean utilisation at swap     : {mean_u:8.3f}",
        f"Mean canister life           : {np.mean(np.diff([0.0] + sim.swap_times))/3600:8.2f} h",
        f"Canisters consumed (sim)     : {sim.canisters_used:8d}",
        f"Canisters to fly (+{spares:.0%} spares): {n_flight:8d}",
        f"LiOH flown                   : {n_flight * c.lioh_mass:8.1f} kg",
        f"Product water released       : {absorbed * water_per_kg_co2():8.2f} kg",
        f"Mean reaction heat load      : {heat_w:8.1f} W",
        f"Peak / mean ppCO2            : {sim.ppco2.max():.3f} / {np.mean(sim.ppco2):.3f} mmHg",
        f"Mass balance error           : "
        f"{abs(sim.co2_generated - sim.co2_absorbed - sim.co2_cabin_delta):.2e} kg",
    ]
    return "\n".join(lines), sim


if __name__ == "__main__":
    text, _ = report(Mission(), Canister())
    print(text)
