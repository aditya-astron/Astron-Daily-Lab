# Astron Built this AP
# Project 019 | RF link budgets | Mars X-band downlink over a synodic cycle
# Author: Aditya (aditya-astron)
# Date: 2026-11-07 IST
"""Checks against published ephemeris events, closed-form coding limits and
hand-computed link terms.

Run: python3 -m pytest -q   (or python3 test_mars_link.py)
"""
import numpy as np

from ephem import AU_KM, date_from_jd, heliocentric, jd_from_date, mars_geometry, solve_kepler
from link import (K_BOLTZ_DB, RATE_LADDER_BPS, XBandLink, biawgn_capacity, ebn0_limit_biawgn,
                  ebn0_limit_shannon, ebn0_uncoded_bpsk, elevation_deg, pass_volume,
                  quantise_rate)


def _event(start, days, key, fn=np.argmin, step=0.01):
    jd = jd_from_date(*start) + np.arange(0.0, days, step)
    g = mars_geometry(jd)
    i = fn(np.abs(g[key]) if key == "dlon_deg" else g[key])
    return jd[i], {k: v[i] for k, v in g.items()}


# --------------------------------------------------------------------------- ephemeris

def test_julian_date_round_trip():
    assert jd_from_date(2000, 1, 1, 12.0) == 2451545.0           # J2000.0
    assert jd_from_date(1957, 10, 4, 19.44) == 2436116.31        # Meeus example 7.a
    for jd in (2451545.0, 2461351.5, 2462000.25):
        y, m, d = date_from_jd(jd)
        assert abs(jd_from_date(y, m, d) - jd) < 1e-8


def test_kepler_solver_residual():
    for e in (0.0, 0.0167, 0.0934, 0.7, 0.97):
        M = np.linspace(-np.pi, np.pi, 721)
        E = solve_kepler(M, e)
        assert np.max(np.abs(E - e * np.sin(E) - M)) < 1e-12


def test_earth_perihelion_aphelion():
    # Earth is closest to the Sun in early January and farthest in early July.
    jd = jd_from_date(2027, 1, 1) + np.arange(0.0, 365.0, 0.05)
    r = np.linalg.norm(heliocentric("earth", jd), axis=-1)
    assert abs(r.min() - 0.98329) < 2e-4 and abs(r.max() - 1.01671) < 2e-4
    y, m, d = date_from_jd(jd[np.argmin(r)])
    assert m == 1 and 1 <= d <= 6        # EMB perihelion lags Earth's by a day or so


def test_mars_2025_close_approach_and_opposition():
    # Published: closest approach 2025-01-12 at 96.08 million km, opposition 2025-01-16.
    jd, g = _event((2024, 12, 15), 60, "range_au")
    assert abs(jd - jd_from_date(2025, 1, 12)) < 1.0
    assert abs(g["range_au"] * AU_KM / 1e6 - 96.08) < 0.15
    jd, _ = _event((2024, 12, 15), 60, "dlon_deg")
    assert abs(jd - jd_from_date(2025, 1, 16, 2.6)) < 1.0


def test_mars_2026_solar_conjunction():
    # Published: Mars solar conjunction 2026-01-09.
    jd, g = _event((2025, 12, 15), 50, "sep_deg")
    assert abs(jd - jd_from_date(2026, 1, 9)) < 1.0
    assert g["sep_deg"] < 1.5


def test_mars_2027_close_approach():
    # Published: opposition 2027-02-19, closest approach 2027-02-20 at 101.4 million km.
    jd, g = _event((2027, 1, 15), 60, "range_au")
    assert abs(jd - jd_from_date(2027, 2, 20)) < 1.0
    assert abs(g["range_au"] * AU_KM / 1e6 - 101.4) < 0.2
    jd, _ = _event((2027, 1, 15), 60, "dlon_deg")
    assert abs(jd - jd_from_date(2027, 2, 19)) < 1.0


def test_mars_range_bounds_and_light_time():
    jd = jd_from_date(2026, 1, 1) + np.arange(0.0, 15 * 365.0, 1.0)
    g = mars_geometry(jd)
    # Physical bounds a_M(1-e_M) - a_E(1+e_E) and a_M(1+e_M) + a_E(1+e_E).
    assert 0.372 < g["range_au"].min() and g["range_au"].max() < 2.684
    assert np.allclose(g["owlt_s"], g["range_au"] * 499.004784, rtol=1e-8)
    # Small-angle sanity: at opposition SEP is near 180, at conjunction near 0.
    assert g["sep_deg"].max() > 170 and g["sep_deg"].min() < 2


# --------------------------------------------------------------------------- coding

def test_uncoded_bpsk_threshold():
    # Textbook values: 9.59 dB at 1e-5, 10.53 dB at 1e-6.
    assert abs(ebn0_uncoded_bpsk(1e-5) - 9.5880) < 1e-3
    assert abs(ebn0_uncoded_bpsk(1e-6) - 10.5301) < 1e-3


def test_bpsk_monte_carlo_matches_theory():
    rng = np.random.default_rng(19)
    ebn0_db = 6.0
    n = 2_000_000
    bits = rng.integers(0, 2, n)
    x = 1.0 - 2.0 * bits
    sigma = np.sqrt(1.0 / (2.0 * 10 ** (ebn0_db / 10)))
    ber = np.mean((x + sigma * rng.standard_normal(n) < 0) != (bits == 0))
    # Inverting the threshold function at the simulated BER must give back ~6 dB.
    assert abs(ebn0_uncoded_bpsk(ber) - ebn0_db) < 0.05


def test_biawgn_capacity_limits():
    assert biawgn_capacity(30.0) > 0.99999
    assert biawgn_capacity(-30.0) < 2e-3
    # Low-SNR slope: C ~ (Es/N0) / ln 2, the same as the Gaussian channel.
    snr = 10 ** (-25 / 10)
    assert abs(biawgn_capacity(-25.0) / (snr / np.log(2)) - 1) < 0.01
    # Binary-input never beats the Gaussian-input capacity.
    for es in np.linspace(-10, 10, 21):
        assert biawgn_capacity(es) <= 0.5 * np.log2(1 + 2 * 10 ** (es / 10)) + 1e-12


def test_coding_limits():
    assert abs(ebn0_limit_biawgn(0.5) - 0.187) < 0.002       # classic rate-1/2 BI-AWGN limit
    assert abs(ebn0_limit_biawgn(1 / 3) - (-0.495)) < 0.005
    assert abs(ebn0_limit_biawgn(1 / 6) - (-1.08)) < 0.02     # CCSDS turbo r=1/6 regime
    assert abs(ebn0_limit_shannon(1.0)) < 1e-12               # rate 1/2, real channel: 0 dB
    assert abs(ebn0_limit_shannon(1e-9) - 10 * np.log10(np.log(2))) < 1e-6   # -1.59 dB
    assert XBandLink().ebn0_req_db >= ebn0_limit_biawgn(0.5)


# --------------------------------------------------------------------------- link

def test_link_terms_by_hand():
    link = XBandLink()
    lam = 299_792_458.0 / 8.42e9
    assert abs(K_BOLTZ_DB - (-228.599)) < 1e-3
    # Free-space loss at exactly 1 AU, 8.42 GHz, computed independently.
    fsl = 20 * np.log10(4 * np.pi * 149_597_870_700.0 / lam)
    assert abs(link.free_space_loss_db(AU_KM) - fsl) < 1e-9
    assert abs(fsl - 274.45) < 0.01
    # 3 m dish, 55 %: G = 0.55 (pi D / lambda)^2.
    assert abs(link.hga_gain_dbi() - 10 * np.log10(0.55 * (np.pi * 3.0 / lam) ** 2)) < 1e-9
    # A 34 m aperture at 74 % reproduces the 68.24 dBi used for the ground station.
    assert abs(10 * np.log10(0.74 * (np.pi * 34.0 / lam) ** 2) - link.gr_dbi) < 0.02
    # Zenith Tsys: 17 + 275 (1 - 10^-0.004) + 2.725 * 10^-0.004
    L = 10 ** 0.004
    assert abs(link.tsys_k(90.0) - (17 + 275 * (1 - 1 / L) + 2.725 / L)) < 1e-9


def test_rate_scales_inverse_square_with_range():
    link = XBandLink()
    r1 = link.max_rate_bps(1.0 * AU_KM, 45.0)
    r2 = link.max_rate_bps(2.0 * AU_KM, 45.0)
    assert abs(r1 / r2 - 4.0) < 1e-9
    # Lower elevation always costs rate (more loss and more noise).
    assert link.max_rate_bps(AU_KM, 10.0) < link.max_rate_bps(AU_KM, 30.0) < r1 * 1.01


def test_budget_identity():
    link = XBandLink()
    d, el = 1.5 * AU_KM, 40.0
    lhs = 10 * np.log10(link.max_rate_bps(d, el))
    rhs = link.pr_n0_dbhz(d, el) - link.ebn0_req_db - link.margin_db
    assert abs(lhs - rhs) < 1e-9


def test_quantise_rate():
    assert quantise_rate(400.0) == 0.0
    assert quantise_rate(500.0) == 500.0
    assert quantise_rate(1_999.0) == 1000.0
    assert quantise_rate(1e9) == RATE_LADDER_BPS[-1]


def test_elevation_geometry():
    lat = 35.426
    assert abs(elevation_deg(lat, 10.0, 0.0) - (90 - abs(lat - 10))) < 1e-9
    # Object on the celestial equator sets at H = 90 deg from any latitude.
    assert abs(elevation_deg(lat, 0.0, 90.0)) < 1e-9


def test_pass_volume_against_fixed_rate_hand_calc():
    link = XBandLink()
    # Dec = +35 at Goldstone: the target is above 10 deg for the whole 8 h track.
    d = 2.0 * AU_KM
    vdr, fixed, hours = pass_volume(link, d, 35.0, 35.426)
    assert abs(hours - 8.0) < 1e-9
    el_edge = elevation_deg(35.426, 35.0, 15.0 * (4.0 - 1.0 / 60.0))   # first midpoint sample
    expect = quantise_rate(link.max_rate_bps(d, el_edge)) * 8 * 3600
    assert abs(fixed - expect) < 1e-6
    assert vdr >= fixed
    # Target far south never rises above 10 deg.
    assert pass_volume(link, d, -60.0, 35.426) == (0.0, 0.0, 0.0)


if __name__ == "__main__":
    import sys
    tests = [(k, v) for k, v in dict(globals()).items() if k.startswith("test_")]
    for name, fn in tests:
        fn()
        print(f"ok  {name}")
    print(f"{len(tests)} passed")
    sys.exit(0)
