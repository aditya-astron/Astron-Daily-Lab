# Astron Built this AP
# Project 019 | RF link budgets | Mars X-band downlink over a synodic cycle
# Author: Aditya (aditya-astron)
# Date: 2026-11-07 IST
"""Run a Mars orbiter's X-band downlink through one full synodic period.

For every day: Earth-Mars geometry, the best supportable rate at the top of
the Goldstone track, and the data volume of one 8 h track with fixed and
variable data rate. Prints the budget at the two range extremes, a monthly
summary, conjunction outage and cycle totals.

Usage: python3 mars_link.py [YYYY-MM-DD start] [days]
"""
import sys

import numpy as np

from ephem import AU_KM, date_from_jd, jd_from_date, mars_geometry
from link import XBandLink, biawgn_capacity, db, ebn0_limit_biawgn, quantise_rate, pass_volume

GOLDSTONE_LAT = 35.426


def fmt_date(jd):
    y, m, d = date_from_jd(jd)
    return f"{y:04d}-{m:02d}-{int(d):02d}"


def fmt_rate(bps):
    if bps >= 1e6:
        return f"{bps / 1e6:.3f} Mbps"
    if bps >= 1e3:
        return f"{bps / 1e3:.1f} kbps"
    return f"{bps:.0f} bps"


def budget_table(link, range_km, el_deg):
    rows = [
        ("Transmit power", db(link.pt_w), "dBW"),
        ("Tx circuit loss", -link.l_tx_db, "dB"),
        (f"HGA gain ({link.hga_diam_m:.1f} m, eff {link.hga_eff})", link.hga_gain_dbi(), "dBi"),
        ("HGA pointing loss", -link.l_pt_sc_db, "dB"),
        ("EIRP", link.eirp_dbw(), "dBW"),
        (f"Free-space loss ({range_km / AU_KM:.3f} AU)", -link.free_space_loss_db(range_km), "dB"),
        (f"Atmosphere ({el_deg:.0f} deg el)", -link.atm_loss_db(el_deg), "dB"),
        ("Polarisation loss", -link.l_pol_db, "dB"),
        ("DSN 34 m gain", link.gr_dbi, "dBi"),
        ("Ground pointing loss", -link.l_pt_gs_db, "dB"),
        (f"Tsys ({link.tsys_k(el_deg):.1f} K)", db(link.tsys_k(el_deg)), "dBK"),
        ("Pr/N0", link.pr_n0_dbhz(range_km, el_deg), "dB-Hz"),
        ("Required Eb/N0 (r=1/2 code)", link.ebn0_req_db, "dB"),
        ("Margin", link.margin_db, "dB"),
        ("Max information rate", db(link.max_rate_bps(range_km, el_deg)), "dB-bps"),
    ]
    for name, val, unit in rows:
        print(f"  {name:<38s} {val:9.2f} {unit}")
    r = link.max_rate_bps(range_km, el_deg)
    print(f"  {'-> rate ladder':<38s} {fmt_rate(quantise_rate(r)):>15s} (unquantised {fmt_rate(r)})")


def run(start=(2026, 11, 7), days=None):
    link = XBandLink()
    jd0 = jd_from_date(*start)
    if days is None:
        # Mean synodic period from the two mean motions in the ephemeris.
        days = int(np.ceil(1.0 / (1.0 / 365.256 - 1.0 / 686.980)))
    jd = jd0 + np.arange(days)
    g = mars_geometry(jd)
    rng_km = g["range_au"] * AU_KM

    peak_el = 90.0 - np.abs(GOLDSTONE_LAT - g["dec_deg"])
    outage = g["sep_deg"] < link.sep_outage_deg
    vdr = np.zeros(days)
    fixed = np.zeros(days)
    peak_rate = np.zeros(days)
    for i in range(days):
        if outage[i]:
            continue
        vdr[i], fixed[i], _ = pass_volume(link, rng_km[i], g["dec_deg"][i], GOLDSTONE_LAT)
        peak_rate[i] = quantise_rate(link.max_rate_bps(rng_km[i], peak_el[i]))

    print(f"Mars X-band downlink, {fmt_date(jd[0])} to {fmt_date(jd[-1])} ({days} days)")
    print(f"Coding: BI-AWGN r=1/2 limit {ebn0_limit_biawgn(0.5):.3f} dB, "
          f"design threshold {link.ebn0_req_db:.2f} dB\n")

    for label, i in (("closest", np.argmin(g["range_au"])), ("farthest", np.argmax(g["range_au"]))):
        print(f"Link budget at {label} approach, {fmt_date(jd[i])}, "
              f"OWLT {g['owlt_s'][i] / 60:.1f} min, peak el {peak_el[i]:.0f} deg:")
        budget_table(link, rng_km[i], peak_el[i])
        print()

    print(" month     range AU     SEP deg   peak rate     Gbit/day VDR   Gbit/day fixed")
    months = np.array([fmt_date(j)[:7] for j in jd])
    for m in dict.fromkeys(months):
        k = months == m
        print(f" {m}  {g['range_au'][k].min():5.3f}-{g['range_au'][k].max():5.3f}"
              f"  {g['sep_deg'][k].min():6.1f}  {fmt_rate(peak_rate[k].max()):>12s}"
              f"  {vdr[k].mean() / 1e9:12.2f}  {fixed[k].mean() / 1e9:14.2f}")

    i_conj = np.argmin(g["sep_deg"])
    o = np.flatnonzero(outage)
    print(f"\nSolar conjunction {fmt_date(jd[i_conj])}, min SEP {g['sep_deg'][i_conj]:.2f} deg")
    if o.size:
        print(f"Outage (SEP < {link.sep_outage_deg:.0f} deg): {fmt_date(jd[o[0]])} to "
              f"{fmt_date(jd[o[-1]])}, {o.size} days")
    print(f"Cycle volume, one 8 h track/day: VDR {vdr.sum() / 8e12:.2f} TB, "
          f"fixed rate {fixed.sum() / 8e12:.2f} TB, VDR gain {100 * (vdr.sum() / fixed.sum() - 1):.1f} %")
    print(f"Best/worst daily volume ratio: {vdr.max() / vdr[vdr > 0].min():.1f}x  "
          f"(range^2 ratio {(g['range_au'].max() / g['range_au'].min()) ** 2:.1f}x)")
    return {"jd": jd, "geom": g, "vdr": vdr, "fixed": fixed, "outage": outage}


if __name__ == "__main__":
    start = tuple(int(x) for x in sys.argv[1].split("-")) if len(sys.argv) > 1 else (2026, 11, 7)
    days = int(sys.argv[2]) if len(sys.argv) > 2 else None
    run(start, days)
