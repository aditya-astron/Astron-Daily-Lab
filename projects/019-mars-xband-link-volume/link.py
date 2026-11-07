# Astron Built this AP
# Project 019 | RF link budgets | Mars X-band downlink over a synodic cycle
# Author: Aditya (aditya-astron)
# Date: 2026-11-07 IST
"""X-band Mars-to-DSN link budget, coding thresholds and pass data volume.

Received power-to-noise density (all in dB):

    Pr/N0 = Pt + Gt - Ltx - Lpt_sc - Lfs - Latm(el) - Lpol + Gr - Lpt_gs - 10 log10(k Tsys(el))
    Lfs   = 20 log10(4 pi d / lambda)
    Latm  = A_zenith / sin(el)                                  (flat-earth slab)
    Tsys  = T_rx + T_phys (1 - 1/Latm) + T_cmb / Latm           (absorbing slab)

The highest supportable information rate is

    Rb_max = Pr/N0 - (Eb/N0)_req - margin    [dB-bps]

rounded down to a power-of-two rate ladder, the way a flight radio's rate
table works.
"""
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq
from scipy.special import erfcinv

K_BOLTZ_DB = 10 * np.log10(1.380649e-23)   # -228.60 dBW/K/Hz
C_M_S = 299_792_458.0
T_CMB = 2.725


def db(x):
    return 10.0 * np.log10(x)


def undb(x):
    return 10.0 ** (np.asarray(x) / 10.0)


# --------------------------------------------------------------------------- coding

def ebn0_uncoded_bpsk(ber):
    """Eb/N0 [dB] for coherent BPSK at a target BER: Pb = Q(sqrt(2 Eb/N0))."""
    return db(erfcinv(2.0 * ber) ** 2)


_GH_X, _GH_W = np.polynomial.hermite.hermgauss(80)


def biawgn_capacity(esn0_db):
    """Capacity [bit/symbol] of binary antipodal input on AWGN.

    With unit amplitude and noise variance s2 = 1/(2 Es/N0),

        C = 1 - E_y[ log2(1 + exp(-2 y / s2)) ],   y ~ N(1, s2)

    The expectation is done with 80-point Gauss-Hermite quadrature.
    """
    s2 = 1.0 / (2.0 * undb(esn0_db))
    y = 1.0 + np.sqrt(2.0 * s2) * _GH_X
    llr = 2.0 * y / s2
    # log2(1 + exp(-llr)) without overflow for large negative llr
    term = np.logaddexp(0.0, -llr) / np.log(2.0)
    return 1.0 - np.sum(_GH_W * term) / np.sqrt(np.pi)


def ebn0_limit_biawgn(rate):
    """Smallest Eb/N0 [dB] at which a rate-R binary code can be error-free."""
    f = lambda ebn0: biawgn_capacity(ebn0 + db(rate)) - rate
    return brentq(f, -2.0, 10.0, xtol=1e-10)


def ebn0_limit_shannon(rate_bps_per_hz):
    """Unconstrained-input limit: Eb/N0 = (2^eta - 1) / eta."""
    eta = rate_bps_per_hz
    return db((2.0 ** eta - 1.0) / eta)


# --------------------------------------------------------------------------- link

@dataclass
class XBandLink:
    freq_hz: float = 8.42e9           # DSN X-band deep-space downlink allocation
    pt_w: float = 100.0               # TWTA RF output
    l_tx_db: float = 1.0              # waveguide, diplexer, switches
    hga_diam_m: float = 3.0
    hga_eff: float = 0.55
    l_pt_sc_db: float = 0.3           # spacecraft HGA pointing
    l_pol_db: float = 0.1
    gr_dbi: float = 68.24             # 34 m BWG X-band gain at the feed (DSN 810-005 class)
    l_pt_gs_db: float = 0.1
    t_rx_k: float = 17.0              # LNA + feed + antenna, excluding atmosphere
    atm_zenith_db: float = 0.04       # clear-sky X-band zenith attenuation, Goldstone
    t_phys_k: float = 275.0           # effective physical temperature of the absorbing layer
    ebn0_req_db: float = 1.19         # BI-AWGN r=1/2 limit (0.19 dB) + 1.0 dB finite-length gap
    margin_db: float = 3.0
    min_el_deg: float = 10.0
    sep_outage_deg: float = 3.0       # solar scintillation kills the link inside this

    @property
    def wavelength_m(self):
        return C_M_S / self.freq_hz

    def hga_gain_dbi(self):
        return db(self.hga_eff * (np.pi * self.hga_diam_m / self.wavelength_m) ** 2)

    def eirp_dbw(self):
        return db(self.pt_w) - self.l_tx_db + self.hga_gain_dbi() - self.l_pt_sc_db

    def free_space_loss_db(self, range_km):
        return 20.0 * np.log10(4.0 * np.pi * np.asarray(range_km) * 1e3 / self.wavelength_m)

    def atm_loss_db(self, el_deg):
        return self.atm_zenith_db / np.sin(np.radians(el_deg))

    def tsys_k(self, el_deg):
        L = undb(self.atm_loss_db(el_deg))
        return self.t_rx_k + self.t_phys_k * (1.0 - 1.0 / L) + T_CMB / L

    def pr_n0_dbhz(self, range_km, el_deg):
        return (self.eirp_dbw() - self.free_space_loss_db(range_km) - self.atm_loss_db(el_deg)
                - self.l_pol_db + self.gr_dbi - self.l_pt_gs_db
                - K_BOLTZ_DB - db(self.tsys_k(el_deg)))

    def max_rate_bps(self, range_km, el_deg):
        """Continuous (unquantised) information rate that closes with margin."""
        return undb(self.pr_n0_dbhz(range_km, el_deg) - self.ebn0_req_db - self.margin_db)


# Rate table: 500 bps ... 8.192 Mbps in powers of two.
RATE_LADDER_BPS = 1000.0 * 2.0 ** np.arange(-1, 14)


def quantise_rate(rate_bps):
    """Highest ladder rate not above rate_bps (0 if none fits)."""
    rate_bps = np.asarray(rate_bps, dtype=float)
    idx = np.searchsorted(RATE_LADDER_BPS, rate_bps, side="right") - 1
    return np.where(idx >= 0, RATE_LADDER_BPS[np.clip(idx, 0, None)], 0.0)


def elevation_deg(lat_deg, dec_deg, hour_angle_deg):
    """sin el = sin(lat) sin(dec) + cos(lat) cos(dec) cos(H)."""
    lat, dec, H = np.radians(lat_deg), np.radians(dec_deg), np.radians(hour_angle_deg)
    return np.degrees(np.arcsin(np.sin(lat) * np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(H)))


def pass_volume(link, range_km, dec_deg, lat_deg, track_hours=8.0, dt_min=2.0):
    """Data volume [bit] of one track centred on transit.

    Returns (vdr_bits, fixed_bits, track_hours_used).

    vdr:   variable data rate, the radio steps along the ladder as elevation
           changes (a time-tagged rate table uploaded ahead of the pass).
    fixed: one rate for the whole track, set by the worst elevation in it.
    """
    n = int(round(track_hours * 60.0 / dt_min))
    t = (np.arange(n) + 0.5) * dt_min / 60.0 - 0.5 * track_hours   # midpoint rule
    el = elevation_deg(lat_deg, dec_deg, 15.0 * t)   # 15 deg/h; Mars' own motion is ~0.03 deg/h
    up = el >= link.min_el_deg
    if not np.any(up):
        return 0.0, 0.0, 0.0
    dt_s = dt_min * 60.0
    rates = quantise_rate(link.max_rate_bps(range_km, np.where(up, el, 90.0)))
    vdr = np.sum(rates[up]) * dt_s
    fixed_rate = quantise_rate(link.max_rate_bps(range_km, el[up].min()))
    fixed = fixed_rate * np.count_nonzero(up) * dt_s
    return float(vdr), float(fixed), np.count_nonzero(up) * dt_min / 60.0
