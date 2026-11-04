# Astron Built this AP
# Project 012 | Radiation shielding | Solar proton event dose behind aluminium vs polyethylene shielding
# Author: Aditya (aditya-astron)
# Date: 2026-11-04 IST
"""Straight-ahead CSDA transport of a solar proton spectrum through a slab.

A proton leaving the slab with energy E_out entered with
E_in = E(R(E_out) + x), x the areal density along its path. Protons are
conserved, so phi_out dE_out = phi_in dE_in, and

    dE_in/dE_out = S(E_in) / S(E_out)

per layer. Dose in a thin tissue (water) layer behind the stack is

    D = C * integral phi_out(E) S_water(E) dE,   C = 1.602e-10 Gy per MeV/g.

For an isotropic omnidirectional fluence F, a point behind an infinite slab
sees F/(4 pi) per steradian through path t/mu from the shielded hemisphere:

    D_iso(t) = 1/2 * integral_0^1 D_normal(t / mu) dmu

the other hemisphere is taken as shadowed by the rest of the vehicle.

Nuclear reactions (primary attenuation, secondary neutrons) are ignored. They
remove a few % of primaries per 10 g/cm^2 and matter mainly at large depths.
"""
import numpy as np

from stopping import E_MAX, E_MIN, WATER, ALUMINIUM, POLYETHYLENE

MEV_PER_G_TO_GY = 1.602176634e-10
E_CUT = 10.0


def king_aug1972(E):
    """King (1974) exponential fit to the August 1972 event.

    Integral fluence Phi(>E) = 7.9e9 exp((30 - E)/26.5) protons/cm^2,
    returned here as the differential fluence (per MeV). The fit is cut at
    10 MeV: extrapolating it lower would let a guess dominate the
    lightly-shielded dose."""
    E = np.asarray(E, dtype=float)
    return np.where(E >= E_CUT, 7.9e9 / 26.5 * np.exp((30.0 - E) / 26.5), 0.0)


def king_integral(E):
    return 7.9e9 * np.exp((30.0 - E) / 26.5)


def transmit(E_out, layers):
    """Map exit energies back through the layers.

    layers: list of (Material, areal density g/cm^2), outermost first.
    Returns (E_in, dE_in/dE_out)."""
    E = np.asarray(E_out, dtype=float)
    jac = np.ones_like(E)
    for mat, x in reversed(layers):
        E_enter = mat.energy(mat.range(E) + x)
        jac *= mat.S(E_enter) / mat.S(E)
        E = E_enter
    return E, jac


def spectrum_out(spectrum, layers, E_out):
    E_in, jac = transmit(E_out, layers)
    phi = spectrum(E_in) * jac
    phi[E_in >= E_MAX * 0.999] = 0.0     # outside the range tables
    return phi


def dose_normal(spectrum, layers, n=1500):
    """Dose (Gy) in thin tissue behind the layers for a normally incident fluence."""
    # The spectrum steps to zero at E_CUT; put grid nodes either side of it so
    # the unshielded case integrates the step exactly.
    E = np.union1d(np.geomspace(E_MIN, E_MAX, n), [E_CUT * (1 - 1e-12), E_CUT])
    f = spectrum_out(spectrum, layers, E) * WATER.S(E) * E   # d/dlnE
    return MEV_PER_G_TO_GY * np.trapezoid(f, np.log(E))


def dose_isotropic(spectrum, material, t, n_mu=48):
    """Dose (Gy) behind a slab of areal density t for an isotropic omnidirectional fluence."""
    if t == 0:
        return 0.5 * dose_normal(spectrum, [])
    x, w = np.polynomial.legendre.leggauss(n_mu)
    mu, w = 0.5 * (x + 1), 0.5 * w
    d = np.array([dose_normal(spectrum, [(material, t / m)]) for m in mu])
    return 0.5 * np.dot(w, d)


def areal_density_for(dose_limit, spectrum, material, lo=0.1, hi=60.0):
    """Bisect on log t for the slab that brings D_iso down to dose_limit."""
    for _ in range(50):
        mid = np.sqrt(lo * hi)
        if dose_isotropic(spectrum, material, mid) > dose_limit:
            lo = mid
        else:
            hi = mid
    return np.sqrt(lo * hi)


def main():
    mats = [ALUMINIUM, POLYETHYLENE, WATER]
    print("King (1974) fit, August 1972 SPE, isotropic, 2-pi slab exposure")
    print("Thin-tissue dose behind slab, Gy (no nuclear attenuation)\n")
    print(f"{'g/cm^2':>7} " + " ".join(f"{m.name:>13}" for m in mats) + "   HDPE/Al")
    D0 = dose_isotropic(king_aug1972, ALUMINIUM, 0.0)
    print(f"{0:7.1f} " + " ".join(f"{D0:13.2f}" for _ in mats) + f"   {1:7.3f}")
    for t in [0.5, 1, 2, 5, 10, 20, 30]:
        d = [dose_isotropic(king_aug1972, m, t) for m in mats]
        print(f"{t:7.1f} " + " ".join(f"{v:13.3f}" for v in d) + f"   {d[1] / d[0]:7.3f}")

    print("\nAreal density to hold this event to the dose, g/cm^2:")
    for lim in [2.0, 1.0, 0.5, 0.25]:
        ta = areal_density_for(lim, king_aug1972, ALUMINIUM)
        tp = areal_density_for(lim, king_aug1972, POLYETHYLENE)
        print(f"  {lim:4.2f} Gy   Al {ta:6.2f} ({ta / 2.70:5.2f} cm)   "
              f"HDPE {tp:6.2f} ({tp / 0.95:5.2f} cm)   mass saving {1 - tp / ta:5.1%}")


if __name__ == "__main__":
    main()
