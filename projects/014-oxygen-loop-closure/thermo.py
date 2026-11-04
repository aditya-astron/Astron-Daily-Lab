# Astron Built this AP
# Project 014 | ECLSS | Oxygen loop closure: PEM electrolysis and Sabatier CO2 reduction
# Author: Aditya (aditya-astron)
# Date: 2026-11-04 IST
"""Ideal-gas thermochemistry from NIST Shomate fits.

    t   = T / 1000
    Cp  = A + B t + C t^2 + D t^3 + E / t^2                         J/mol/K
    H   = A t + B t^2/2 + C t^3/3 + D t^4/4 - E/t + F               kJ/mol
    S   = A ln t + B t + C t^2/2 + D t^3/3 - E/(2 t^2) + G          J/mol/K

H here is the absolute enthalpy on the standard-formation scale, so
H(298.15) is the enthalpy of formation. G = H - T S. Coefficient sets are
the NIST WebBook ranges that cover 298-1000 K (water vapour's fit starts at
500 K but stays within 0.05 J/mol/K of JANAF down to 298 K).
"""
import math

R = 8.314462618          # J/mol/K
FARADAY = 96485.33212    # C/mol
T_REF = 298.15

#            A           B           C           D           E          F           G          H(f,298)
SHOMATE = {
    "H2":   (33.066178, -11.363417, 11.432816, -2.772874, -0.158558, -9.980797, 172.707974, 0.0),
    "O2":   (31.32234,  -20.23531,  57.86644,  -36.50624, -0.007374, -8.903471, 246.7945,   0.0),
    "CO2":  (24.99735,   55.18696, -33.69137,   7.948387, -0.136638, -403.6075, 228.2431, -393.5224),
    "CH4":  (-0.703029, 108.4773,  -42.52157,   5.862788,  0.678565, -76.84376, 158.7163,  -74.87310),
    "H2O":  (30.09200,   6.832514,  6.793435,  -2.534480,  0.082139, -250.8810, 223.3967, -241.8264),
    "H2O(l)": (-203.6060, 1523.290, -3196.413, 2474.455,   3.855326, -256.5478, -488.7163, -285.8304),
}


def cp(sp, T):
    A, B, C, D, E, *_ = SHOMATE[sp]
    t = T / 1000.0
    return A + B * t + C * t**2 + D * t**3 + E / t**2


def enthalpy(sp, T):
    """Absolute molar enthalpy, J/mol (formation scale)."""
    A, B, C, D, E, F, _, _ = SHOMATE[sp]
    t = T / 1000.0
    return 1000.0 * (A * t + B * t**2 / 2 + C * t**3 / 3 + D * t**4 / 4 - E / t + F)


def entropy(sp, T):
    A, B, C, D, E, _, G, _ = SHOMATE[sp]
    t = T / 1000.0
    return A * math.log(t) + B * t + C * t**2 / 2 + D * t**3 / 3 - E / (2 * t**2) + G


def reaction(nu, T):
    """(dH, dS, dG) of a reaction given stoichiometry {species: nu}, nu < 0 for reactants."""
    dH = sum(n * enthalpy(s, T) for s, n in nu.items())
    dS = sum(n * entropy(s, T) for s, n in nu.items())
    return dH, dS, dH - T * dS


# CO2 + 4 H2 -> CH4 + 2 H2O(g)
SABATIER = {"CO2": -1, "H2": -4, "CH4": 1, "H2O": 2}
# H2O(l) -> H2 + 1/2 O2
ELECTROLYSIS = {"H2O(l)": -1, "H2": 1, "O2": 0.5}


def ln_K(nu, T):
    """Natural log of the standard-state (1 bar) equilibrium constant."""
    return -reaction(nu, T)[2] / (R * T)
