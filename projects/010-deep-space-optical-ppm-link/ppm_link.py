# Astron Built this AP
# Project 010 | Laser/optical communication and link budgets | Photon-counting PPM link for a deep-space optical downlink
# Author: Aditya (aditya-astron)
# Date: 2026-11-03 IST
"""Photon-counting M-PPM downlink from a Mars-bound spacecraft to a 5 m ground telescope.

Three pieces:
  1. Link budget: transmit power -> mean detected signal photons per second.
  2. Poisson PPM channel: exact uncoded symbol error rate and the channel
     capacity (bits per PPM symbol) with background light.
  3. Rate optimisation: for each range, pick the PPM order M that maximises
     throughput C(M) / (M * Ts) with an implementation margin from capacity.

Detector dead time, blocking and timing jitter are ignored; slots are assumed
perfectly synchronised.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.special import logsumexp
from scipy.stats import poisson

H_PLANCK = 6.62607015e-34  # J s
C_LIGHT = 299_792_458.0  # m/s
AU = 1.495978707e11  # m


def db(x: float) -> float:
    return 10.0 * math.log10(x)


def from_db(x_db: float) -> float:
    return 10.0 ** (x_db / 10.0)


# ---------------------------------------------------------------------------
# 1. Link budget
# ---------------------------------------------------------------------------
@dataclass
class OpticalLink:
    wavelength: float = 1550e-9  # m
    p_tx: float = 4.0  # W average laser power
    d_tx: float = 0.22  # m flight transceiver aperture
    d_rx: float = 5.0  # m ground telescope aperture
    rx_obscuration: float = 0.36  # m^2/m^2 area lost to secondary mirror
    tx_optics_db: float = -2.0
    pointing_db: float = -3.0  # mean mispointing loss
    atmosphere_db: float = -1.5  # clear-sky transmission at 30 deg elevation
    rx_optics_db: float = -3.0  # telescope train, filter, fibre coupling
    detector_efficiency: float = 0.5  # SNSPD system detection efficiency
    # Background: sky radiance seen in the receiver field of view.
    sky_radiance: float = 0.1  # W / (m^2 sr um), twilight sky near 1.55 um
    fov_half_angle: float = 50e-6  # rad
    filter_width_um: float = 1e-3  # 1 nm optical filter

    @property
    def photon_energy(self) -> float:
        return H_PLANCK * C_LIGHT / self.wavelength

    @property
    def tx_gain(self) -> float:
        # Diffraction-limited circular aperture gain, (pi D / lambda)^2.
        return (math.pi * self.d_tx / self.wavelength) ** 2

    @property
    def rx_area(self) -> float:
        return math.pi * self.d_rx**2 / 4.0 * (1.0 - self.rx_obscuration)

    @property
    def rx_gain(self) -> float:
        # G = 4 pi A / lambda^2 using the effective (obscured) area.
        return 4.0 * math.pi * self.rx_area / self.wavelength**2

    def free_space_loss(self, r: float) -> float:
        return (self.wavelength / (4.0 * math.pi * r)) ** 2

    def losses(self) -> float:
        return from_db(self.tx_optics_db + self.pointing_db + self.atmosphere_db + self.rx_optics_db)

    def received_power(self, r: float) -> float:
        """Optical power at the detector, W (Friis form)."""
        return self.p_tx * self.tx_gain * self.rx_gain * self.free_space_loss(r) * self.losses()

    def signal_photon_rate(self, r: float) -> float:
        """Detected signal photons per second."""
        return self.received_power(r) * self.detector_efficiency / self.photon_energy

    def background_photon_rate(self) -> float:
        """Detected background photons per second: L * A * Omega * d_lambda."""
        omega = math.pi * self.fov_half_angle**2
        # Radiance is specified at the ground, so only receiver optics attenuate it.
        p_bg = self.sky_radiance * self.rx_area * omega * self.filter_width_um * from_db(self.rx_optics_db)
        return p_bg * self.detector_efficiency / self.photon_energy

    def budget_table(self, r: float) -> list[tuple[str, float]]:
        return [
            ("Transmit power [dBW]", db(self.p_tx)),
            ("Transmit gain [dBi]", db(self.tx_gain)),
            ("Free-space loss [dB]", db(self.free_space_loss(r))),
            ("Receive gain [dBi]", db(self.rx_gain)),
            ("Optics/pointing/atmosphere [dB]", db(self.losses())),
            ("Received power [dBW]", db(self.received_power(r))),
        ]


# ---------------------------------------------------------------------------
# 2. Poisson PPM channel
# ---------------------------------------------------------------------------
def ppm_ser(M: int, ks: float, kb: float, kmax: int | None = None) -> float:
    """Exact symbol error rate of uncoded M-PPM, max-count detection, random tie break.

    Signal slot counts ~ Pois(ks + kb), each of the M-1 empty slots ~ Pois(kb).
    Conditioning on the signal-slot count k and on j empty slots tying at k:

        Pc = sum_k p1(k) sum_j C(M-1, j) pb(k)^j Fb(k-1)^(M-1-j) / (j+1)
           = sum_k p1(k) [Fb(k)^M - Fb(k-1)^M] / (M pb(k))

    where the second line is the binomial identity
    sum_j C(n,j) a^j b^(n-j) / (j+1) = ((a+b)^(n+1) - b^(n+1)) / ((n+1) a).
    """
    if kb == 0.0:
        # Only an erasure (zero signal counts) can cause an error.
        return (1.0 - 1.0 / M) * math.exp(-ks)
    if kmax is None:
        kmax = int(ks + kb + 12.0 * math.sqrt(ks + kb) + 30)
    k = np.arange(kmax + 1)
    p1 = poisson.pmf(k, ks + kb)
    pb = poisson.pmf(k, kb)
    Fb = poisson.cdf(k, kb)
    with np.errstate(divide="ignore", invalid="ignore"):
        x = pb / Fb
        # Fb(k)^M - Fb(k-1)^M = Fb(k)^M * (1 - (1 - x)^M), computed without cancellation.
        diff = Fb**M * -np.expm1(M * np.log1p(-x))
        ratio = np.where(pb > 1e-300, diff / (M * pb), Fb ** (M - 1))
    pc = float(np.sum(p1 * ratio))
    return 1.0 - pc


def ppm_capacity(M: int, ks: float, kb: float, n_mc: int = 200_000, seed: int = 1) -> float:
    """Capacity of the Poisson M-PPM channel in bits per PPM symbol.

    With uniform inputs the slot likelihood ratio is (1 + ks/kb)^k, so

        C = log2 M - E[ log2 sum_j exp(a (k_j - k_0)) ],   a = ln(1 + ks/kb)

    where slot 0 carries the pulse. The expectation is estimated by Monte Carlo.
    With kb = 0 the channel is an M-ary erasure channel: C = (1 - e^-ks) log2 M.
    """
    if kb == 0.0:
        return (1.0 - math.exp(-ks)) * math.log2(M)
    rng = np.random.default_rng(seed)
    counts = rng.poisson(kb, size=(n_mc, M)).astype(float)
    counts[:, 0] += rng.poisson(ks, size=n_mc)
    a = math.log1p(ks / kb)
    llr = a * (counts - counts[:, :1])
    penalty = logsumexp(llr, axis=1) / math.log(2.0)
    return math.log2(M) - float(np.mean(penalty))


def binary_ppm_capacity_exact(ks: float, kb: float, kmax: int | None = None) -> float:
    """Exact mutual information of 2-PPM by summing over all count pairs (k0, k1)."""
    if kmax is None:
        kmax = int(ks + kb + 12.0 * math.sqrt(ks + kb + 1.0) + 30)
    k = np.arange(kmax + 1)
    ps = poisson.pmf(k, ks + kb)
    pb = poisson.pmf(k, kb)
    y_given_0 = np.outer(ps, pb)  # pulse in slot 0
    y_given_1 = np.outer(pb, ps)  # pulse in slot 1
    py = 0.5 * (y_given_0 + y_given_1)
    mi = 0.0
    for pyx in (y_given_0, y_given_1):
        mask = pyx > 0
        mi += 0.5 * np.sum(pyx[mask] * np.log2(pyx[mask] / py[mask]))
    return float(mi)


# ---------------------------------------------------------------------------
# 3. Throughput optimisation over PPM order
# ---------------------------------------------------------------------------
def best_ppm_order(
    signal_rate: float,
    background_rate: float,
    slot: float,
    orders: tuple[int, ...] = (4, 8, 16, 32, 64, 128, 256),
    margin_db: float = 1.0,
    n_mc: int = 40_000,
) -> tuple[int, float, float, float]:
    """Return (M, Ks, Kb, throughput bit/s) maximising C(M) / (M Ts).

    Average power is fixed, so a symbol of M slots collects Ks = rate * M * Ts
    signal photons in its single pulsed slot. margin_db de-rates Ks to stand in
    for the gap between capacity and a practical SCPPM code.
    """
    best = (0, 0.0, 0.0, -1.0)
    kb = background_rate * slot
    for M in orders:
        ks = signal_rate * M * slot * from_db(-margin_db)
        rate = ppm_capacity(M, ks, kb, n_mc=n_mc) / (M * slot)
        if rate > best[3]:
            best = (M, ks, kb, rate)
    return best


def main() -> None:
    link = OpticalLink()
    slot = 0.5e-9

    r_ref = 0.2 * AU
    print(f"Link budget at {r_ref / AU:.1f} AU, lambda = {link.wavelength * 1e9:.0f} nm")
    for name, value in link.budget_table(r_ref):
        print(f"  {name:<34s}{value:9.2f}")
    ns = link.signal_photon_rate(r_ref)
    nb = link.background_photon_rate()
    print(f"  {'Detected signal photons [1/s]':<34s}{ns:9.3e}")
    print(f"  {'Detected background photons [1/s]':<34s}{nb:9.3e}")
    print()

    print("Uncoded 64-PPM symbol error rate, Kb = 0.01 per slot")
    for ks in (2.0, 5.0, 10.0, 15.0):
        print(f"  Ks = {ks:5.1f}   SER = {ppm_ser(64, ks, 0.01):.3e}")
    print()

    print("Best PPM order per range (Ts = 0.5 ns, 1 dB margin to capacity)")
    print("  range[AU]  M     Ks      Kb      bits/photon  throughput")
    for r_au in (0.2, 0.5, 1.0, 1.5, 2.0, 2.7):
        rate_s = link.signal_photon_rate(r_au * AU)
        M, ks, kb, rate = best_ppm_order(rate_s, nb, slot)
        bpp = rate / rate_s
        print(f"  {r_au:8.1f}  {M:3d}  {ks:7.3f}  {kb:6.4f}  {bpp:10.2f}  {rate / 1e6:8.2f} Mb/s")


if __name__ == "__main__":
    main()
