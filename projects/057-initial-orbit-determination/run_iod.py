# Astron Built this AP
# Project 057 | Orbital mechanics / Navigation | Initial orbit determination: Gibbs, Herrick-Gibbs and Gauss angles-only
# Author: Aditya (aditya-astron)
# Date: 2026-11-22 IST
"""IOD study driver: `python3 run_iod.py` prints every table in the README.

1. Textbook checks (Curtis 5.1 Gibbs, Vallado 7-4 Herrick-Gibbs, Curtis
   5.11/5.12 Gauss) and how much of the Gauss example's residual is just
   rounding in its published inputs.
2. Gibbs vs Herrick-Gibbs on a 700 km SSO orbit: truncation error, 10 m
   position noise, and the separation angle where the two swap places.
3. Gauss angles-only on a LEO pass seen by a 1.68 km telescope site at
   24.65 N: exact data, then 1 arcsec Monte Carlo against the sighting
   spacing.
4. The same method on a 3 deg inclined geosynchronous object over hours
   of arc, where Curtis's fixed-point refinement stops converging and a
   six-state Newton correction on the sightlines does not.
"""
from __future__ import annotations

import numpy as np

import iod

D = np.deg2rad
RNG = np.random.default_rng(57)


def textbook_checks():
    print("1. Reference examples")
    r1 = np.array([-294.32, 4265.1, 5986.7])
    r2 = np.array([-1365.5, 3637.6, 6346.8])
    r3 = np.array([-2940.3, 2473.7, 6555.8])
    v = iod.gibbs(r1, r2, r3, 398600.0)
    print(f"   Gibbs (Curtis 5.1)        v2 = {np.round(v, 4)}  book [-6.2174 -4.0122  1.599 ]"
          f"  coplanarity {np.degrees(iod.coplanarity(r1, r2, r3)) * 3600:.2f} arcsec")
    r1 = np.array([3419.85564, 6019.82602, 2784.60022])
    r2 = np.array([2935.91195, 6326.18324, 2660.59584])
    r3 = np.array([2434.95205, 6597.38674, 2521.52311])
    v = iod.herrick_gibbs(r1, r2, r3, 0.0, 76.48, 153.04)
    print(f"   Herrick-Gibbs (Vallado 7-4) v2 = {np.round(v, 6)}  book [-6.441557  3.777625 -1.720642]")

    th, ra, de, t = curtis_gauss_inputs()
    out = gauss_curtis(th, ra, de, t)
    print(f"   Gauss (Curtis 5.11) root |r2| = {out['roots'][0]:.1f} km   r2 = {np.round(out['r2_initial'], 1)}")
    print(f"   Gauss refined (5.12)  {out['iterations']:2d} it  r2 = {np.round(out['r2'], 1)}"
          f"  book [5662.1 6538.  3269. ]")
    print(f"                                v2 = {np.round(out['v2'], 4)}  book [-3.8856  5.1214 -2.2433]")
    spread = []
    for _ in range(200):
        o = gauss_curtis(th + RNG.uniform(-5e-4, 5e-4, 3), ra + RNG.uniform(-5e-5, 5e-5, 3),
                         de + RNG.uniform(-5e-5, 5e-5, 3), t)
        spread.append(np.linalg.norm(o["r2"] - out["r2"]))
    print(f"   rounding the example's inputs by half a printed digit moves r2 by"
          f" {np.median(spread):.2f} km median, {np.max(spread):.2f} km max:"
          f" the ~1 km gap to the book is inside that\n")


def curtis_gauss_inputs():
    th = np.array([44.506, 45.000, 45.499])
    ra = np.array([43.5365, 54.4196, 64.3178])
    de = np.array([-8.78334, -12.0739, -15.1054])
    t = [0.0, 118.104, 237.577]
    return th, ra, de, t


def gauss_curtis(th, ra, de, t):
    R = [iod.site_eci(D(40.0), 1.0, D(x), 6378.0, 0.003353) for x in th]
    L = [iod.radec_to_unit(D(a), D(b)) for a, b in zip(ra, de)]
    return iod.gauss(L, R, t, 398600.0)


# ---------------------------------------------------------------------------

def gibbs_vs_hg(sigma_km=0.010, trials=2000):
    a, e, inc = 7078.137, 0.01, D(98.2)
    r0, v0 = iod.coe_to_rv(a, e, inc, D(40), D(60), D(0))
    period = 2 * np.pi * np.sqrt(a ** 3 / iod.MU)
    print("2. Gibbs vs Herrick-Gibbs, 700 km SSO, equal spacing, sigma = "
          f"{sigma_km * 1e3:.0f} m per axis ({trials} trials)")
    print("   sep    dt     HG trunc   Gibbs exact   Gibbs RMS   HG RMS     (m/s)")
    rows = []
    for sep in [0.25, 0.5, 1, 2, 3, 5, 8, 12, 20, 30, 45]:
        dt = sep / 360 * period
        ts = np.array([0.0, dt, 2 * dt])
        rs = [iod.propagate(r0, v0, t)[0] for t in ts]
        v2 = iod.propagate(r0, v0, dt)[1]
        e_hg = np.linalg.norm(iod.herrick_gibbs(*rs, *ts) - v2)
        e_g0 = np.linalg.norm(iod.gibbs(*rs) - v2)
        eg, eh = [], []
        for _ in range(trials):
            rn = [r + RNG.normal(0, sigma_km, 3) for r in rs]
            eg.append(np.sum((iod.gibbs(*rn) - v2) ** 2))
            eh.append(np.sum((iod.herrick_gibbs(*rn, *ts) - v2) ** 2))
        g, h = np.sqrt(np.mean(eg)), np.sqrt(np.mean(eh))
        rows.append((sep, dt, e_hg, e_g0, g, h))
        print(f"   {sep:5.2f}  {dt:6.1f}  {e_hg * 1e3:9.2e}  {e_g0 * 1e3:9.2e}   {g * 1e3:9.3f}  {h * 1e3:9.3f}")
    # Noise floor for any 3-point velocity estimate of the middle point:
    # central difference (r3 - r1)/(2 dt) has sigma*sqrt(3/2)/dt total RMS.
    sep, dt = rows[2][0], rows[2][1]
    print(f"   central-difference noise floor at {sep} deg: "
          f"{sigma_km * np.sqrt(1.5) / dt * 1e3:.3f} m/s (3-axis RMS)")
    cross = None
    for (s0, _, _, _, g0, h0), (s1, _, _, _, g1, h1) in zip(rows, rows[1:]):
        if h0 <= g0 and h1 > g1:
            # log-log interpolation of the ratio crossing 1
            r0_, r1_ = np.log(h0 / g0), np.log(h1 / g1)
            cross = np.exp(np.log(s0) + (np.log(s1) - np.log(s0)) * (-r0_) / (r1_ - r0_))
    if cross is not None:
        print(f"   Herrick-Gibbs is better below ~{cross:.0f} deg separation, Gibbs above\n")
    else:
        print("   no crossover in the scanned range\n")
    return rows, cross


# ---------------------------------------------------------------------------

SITE_LAT, SITE_ALT, SITE_LON = D(24.65), 1.68, D(72.78)   # high-altitude telescope site


def site_at(t, theta_g0=0.0):
    return iod.site_eci(SITE_LAT, SITE_ALT, theta_g0 + SITE_LON + iod.OMEGA_E * t)


def find_pass(r0, v0, min_peak_deg=50.0, step=10.0, horizon=4 * 86400.0, mask=D(15)):
    """First pass above the mask whose peak elevation exceeds min_peak_deg."""
    ts = np.arange(0.0, horizon, step)
    el = np.array([iod.elevation(iod.propagate(r0, v0, t)[0], site_at(t)) for t in ts])
    above = el > mask
    k = 0
    while k < len(ts):
        if above[k]:
            j = k
            while j < len(ts) and above[j]:
                j += 1
            if np.max(el[k:j]) > D(min_peak_deg):
                kp = k + int(np.argmax(el[k:j]))
                return ts[k], ts[kp], ts[j - 1], np.degrees(el[kp])
            k = j
        k += 1
    raise RuntimeError("no qualifying pass")


def observe(r0, v0, ts, sigma_rad=0.0):
    """Unit sightlines with isotropic angular noise (sigma per sky axis)."""
    L, R = [], []
    for t in ts:
        r = iod.propagate(r0, v0, t)[0]
        Rk = site_at(t)
        u = (r - Rk) / np.linalg.norm(r - Rk)
        if sigma_rad > 0:
            e1 = np.cross(u, [0.0, 0.0, 1.0])
            e1 /= np.linalg.norm(e1)
            e2 = np.cross(u, e1)
            u = u + sigma_rad * (RNG.normal() * e1 + RNG.normal() * e2)
            u /= np.linalg.norm(u)
        L.append(u)
        R.append(Rk)
    return L, R


def gauss_study(r0, v0, t_mid, spacings, sigma_arcsec, trials, label):
    sig = np.radians(sigma_arcsec / 3600.0)
    r2t, v2t = iod.propagate(r0, v0, t_mid)
    a_t = iod.rv_to_coe(r2t, v2t)[0]
    print("   exact sightings: series error, Curtis fixed point, Newton correction;"
          f" then {sigma_arcsec:g} arcsec noise, {trials} trials (Newton), medians")
    print("   spacing  D0      series   fixed point        Newton        |   |dr2|      |dv2|     |da|    fail")
    rows = []
    for sp in spacings:
        ts = [t_mid - sp, t_mid, t_mid + sp]
        L, R = observe(r0, v0, ts)
        d0 = abs(L[0] @ np.cross(L[1], L[2]))
        o = iod.gauss(L, R, ts, refine=None)
        e_ser = np.linalg.norm(o["r2_initial"] - r2t)
        try:
            of = iod.gauss(L, R, ts, refine="fixed-point")
            fp = f"{of['iterations']:3d} it {np.linalg.norm(of['r2'] - r2t) * 1e3:7.1e} m"
            fp_it = of["iterations"]
        except RuntimeError:
            fp, fp_it = "  diverges       ", None
        oc = iod.gauss(L, R, ts, refine="correct")
        e_c = np.linalg.norm(oc["r2"] - r2t)
        dr, dv, da, fails = [], [], [], 0
        for _ in range(trials):
            Ln, Rn = observe(r0, v0, ts, sig)
            try:
                on = iod.gauss(Ln, Rn, ts, refine="correct")
                an = iod.rv_to_coe(on["r2"], on["v2"])[0]
            except (ValueError, RuntimeError, np.linalg.LinAlgError):
                fails += 1
                continue
            if not np.isfinite(an) or an < 0:      # hyperbolic fit: unusable
                fails += 1
                continue
            dr.append(np.linalg.norm(on["r2"] - r2t))
            dv.append(np.linalg.norm(on["v2"] - v2t))
            da.append(abs(an - a_t))
        row = dict(sp=sp, d0=d0, series=e_ser, fp_it=fp_it, newton=e_c, newton_it=oc["iterations"],
                   dr=np.median(dr), dv=np.median(dv), da=np.median(da), fail=fails / trials)
        rows.append(row)
        print(f"   {label(sp):>7}  {d0:6.1e} {e_ser:8.3f} km  {fp}  {oc['iterations']:2d} it {e_c * 1e3:7.1e} m"
              f"  | {row['dr']:8.3f} km {row['dv'] * 1e3:7.2f} m/s {row['da']:8.2f} km {row['fail'] * 100:4.1f}%")
    return rows


def gauss_leo():
    a, e, inc = 6378.137 + 800.0, 0.001, D(65.0)
    r0, v0 = iod.coe_to_rv(a, e, inc, D(100), D(30), D(0))
    t0, tp, t1, elmax = find_pass(r0, v0)
    print(f"3. Gauss angles-only, 800 km x 65 deg LEO from a 24.65 N, 1.68 km site")
    print(f"   pass above 15 deg: {t1 - t0:.0f} s long, peak elevation {elmax:.1f} deg;"
          f" sightings symmetric about culmination")
    half = (t1 - t0) / 2
    spacings = [s for s in [5, 10, 20, 40, 80, 120, 160] if s < min(tp - t0, t1 - tp)]
    rows = gauss_study(r0, v0, tp, spacings, 1.0, 500, lambda s: f"{s:.0f} s")
    print(f"   (half-pass length {half:.0f} s limits the widest spacing)\n")
    return rows


def gauss_geo():
    a = 42164.17
    r0, v0 = iod.coe_to_rv(a, 0.0005, D(3.0), D(80), D(0), D(0))
    print("4. Gauss angles-only, geosynchronous object, i = 3 deg, same site")
    t_mid = 6 * 3600.0
    r2 = iod.propagate(r0, v0, t_mid)[0]
    el = np.degrees(iod.elevation(r2, site_at(t_mid)))
    print(f"   elevation at the middle sighting {el:.1f} deg (a GEO object barely moves on the sky,"
          " so the whole night is usable)")
    rows = gauss_study(r0, v0, t_mid, [600, 1800, 3600, 5400, 7200, 10800, 14400], 1.0, 300,
                       lambda s: f"{s / 60:.0f} min")
    ts = [0.0, t_mid, 2 * t_mid]
    L, R = observe(r0, v0, ts)
    try:
        iod.gauss(L, R, ts, refine=None)
        msg = "solved"
    except ValueError as exc:
        msg = str(exc)
    print(f"   at 6 h spacing (180 deg of arc) the series start fails: {msg}\n")
    return rows


def main():
    np.set_printoptions(suppress=True)
    textbook_checks()
    gibbs_vs_hg()
    gauss_leo()
    gauss_geo()


if __name__ == "__main__":
    main()
