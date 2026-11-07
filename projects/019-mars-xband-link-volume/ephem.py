# Astron Built this AP
# Project 019 | RF link budgets | Mars X-band downlink over a synodic cycle
# Author: Aditya (aditya-astron)
# Date: 2026-11-07 IST
"""Low-precision Earth and Mars ephemeris and the geometry a link budget needs.

Mean Keplerian elements and their secular rates are the JPL approximate
planetary positions (Standish, table 1, valid 1800-2050, J2000 ecliptic).
Mars comes out good to a few arcminutes, which is ~0.001 AU in range: far
below anything that matters to a link budget, but enough to land oppositions
and conjunctions on the right day.

Earth here is the Earth-Moon barycentre. The ~4700 km offset is 3e-5 AU.
"""
import numpy as np

AU_KM = 149_597_870.7
C_KM_S = 299_792.458
J2000_JD = 2_451_545.0
OBLIQUITY_J2000 = np.radians(23.43928)

# a [AU], e, I, L, long. perihelion, long. node [deg]; then rates per Julian century
_ELEMENTS = {
    "earth": ((1.00000261, 0.01671123, -0.00001531, 100.46457166, 102.93768193, 0.0),
              (0.00000562, -0.00004392, -0.01294668, 35999.37244981, 0.32327364, 0.0)),
    "mars": ((1.52371034, 0.09339410, 1.84969142, -4.55343205, -23.94362959, 49.55953891),
             (0.00001847, 0.00007882, -0.00813131, 19140.30268499, 0.44441088, -0.29257343)),
}


def jd_from_date(year, month, day, hour=0.0):
    """Julian date of a Gregorian calendar date (Meeus ch. 7)."""
    if month <= 2:
        year -= 1
        month += 12
    a = year // 100
    b = 2 - a + a // 4
    return (int(365.25 * (year + 4716)) + int(30.6001 * (month + 1))
            + day + hour / 24.0 + b - 1524.5)


def date_from_jd(jd):
    """Inverse of jd_from_date, returns (year, month, day) with fractional day."""
    jd += 0.5
    z = int(jd)
    f = jd - z
    alpha = int((z - 1867216.25) / 36524.25)
    a = z + 1 + alpha - alpha // 4
    b = a + 1524
    c = int((b - 122.1) / 365.25)
    d = int(365.25 * c)
    e = int((b - d) / 30.6001)
    day = b - d - int(30.6001 * e) + f
    month = e - 1 if e < 14 else e - 13
    year = c - 4716 if month > 2 else c - 4715
    return year, month, day


def solve_kepler(M, e, tol=1e-14):
    """Eccentric anomaly from mean anomaly, Newton from E0 = M + e sin M."""
    M = np.asarray(M, dtype=float)
    E = M + e * np.sin(M)
    for _ in range(50):
        dE = (E - e * np.sin(E) - M) / (1.0 - e * np.cos(E))
        E = E - dE
        if np.all(np.abs(dE) < tol):
            break
    return E


def heliocentric(body, jd):
    """Heliocentric J2000 ecliptic position [AU] of `body` at Julian date(s) jd.

    Returns an array of shape (..., 3).
    """
    el0, rate = _ELEMENTS[body]
    T = (np.asarray(jd, dtype=float) - J2000_JD) / 36525.0
    a, e, inc, L, varpi, node = (x0 + dx * T for x0, dx in zip(el0, rate))
    inc, L, varpi, node = np.radians(inc), np.radians(L), np.radians(varpi), np.radians(node)
    argp = varpi - node
    M = np.remainder(L - varpi + np.pi, 2 * np.pi) - np.pi
    E = solve_kepler(M, e)

    # Perifocal coordinates, then rotate by node, inclination and argument of perihelion.
    xp = a * (np.cos(E) - e)
    yp = a * np.sqrt(1.0 - e * e) * np.sin(E)
    cw, sw = np.cos(argp), np.sin(argp)
    cO, sO = np.cos(node), np.sin(node)
    ci, si = np.cos(inc), np.sin(inc)
    x = (cw * cO - sw * sO * ci) * xp + (-sw * cO - cw * sO * ci) * yp
    y = (cw * sO + sw * cO * ci) * xp + (-sw * sO + cw * cO * ci) * yp
    z = (sw * si) * xp + (cw * si) * yp
    return np.stack([x, y, z], axis=-1)


def mars_geometry(jd):
    """Geometry of Mars as seen from Earth.

    Returns dict with range [AU], Sun-Earth-Probe angle [deg], one-way light
    time [s], geocentric declination [deg], and the heliocentric longitude
    difference Mars - Earth [deg] (zero at opposition).
    """
    rE = heliocentric("earth", jd)
    rM = heliocentric("mars", jd)
    rho = rM - rE
    rng = np.linalg.norm(rho, axis=-1)

    # SEP: angle at Earth between the Sun direction (-rE) and Mars.
    cos_sep = np.sum(-rE * rho, axis=-1) / (np.linalg.norm(rE, axis=-1) * rng)
    sep = np.degrees(np.arccos(np.clip(cos_sep, -1.0, 1.0)))

    # Ecliptic -> equatorial for the declination the ground station sees.
    ce, se = np.cos(OBLIQUITY_J2000), np.sin(OBLIQUITY_J2000)
    z_eq = se * rho[..., 1] + ce * rho[..., 2]
    dec = np.degrees(np.arcsin(z_eq / rng))

    dlon = np.degrees(np.arctan2(rM[..., 1], rM[..., 0]) - np.arctan2(rE[..., 1], rE[..., 0]))
    dlon = (dlon + 180.0) % 360.0 - 180.0

    return {"range_au": rng, "sep_deg": sep, "owlt_s": rng * AU_KM / C_KM_S,
            "dec_deg": dec, "dlon_deg": dlon}
