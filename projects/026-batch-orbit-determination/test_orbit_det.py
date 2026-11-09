# Astron Built this AP
# Project 026 | Navigation | Batch least-squares orbit determination from ground-station range and range-rate
# Author: Aditya (aditya-astron)
# Date: 2026-11-09 IST
"""Verification of the batch orbit determination. Run: python3 test_orbit_det.py (or pytest)."""
import numpy as np
from scipy.stats import chi2

import orbit_det as od

OBS = od.schedule(od.X_TRUE, od.ARC_S, od.DT_S)


def test_ellipsoid_points():
    assert np.allclose(od.geodetic_to_ecef(0, 0, 0), [od.RE, 0, 0], atol=1e-9)
    b = od.RE * (1 - od.F_WGS)                      # polar radius 6356.752 km
    assert np.allclose(od.geodetic_to_ecef(90, 0, 0), [0, 0, b], atol=1e-9)
    assert abs(b - 6356.752314) < 1e-6


def test_j2_equatorial_gravity():
    # In the equatorial plane |a| = mu/r^2 (1 + 1.5 J2 (Re/r)^2), purely radial.
    r = od.RE + 500.0
    a = od.accel(np.array([r, 0.0, 0.0]))
    ref = od.MU / r**2 * (1 + 1.5 * od.J2 * (od.RE / r) ** 2)
    assert abs(-a[0] - ref) < 1e-15 and abs(a[1]) < 1e-18 and abs(a[2]) < 1e-18


def test_gravity_gradient_vs_finite_difference():
    rng = np.random.default_rng(1)
    for _ in range(5):
        r = rng.normal(size=3)
        r *= (od.RE + 300 + 1000 * rng.random()) / np.linalg.norm(r)
        G = od.gravity_gradient(r)
        h = 1e-3
        Gfd = np.column_stack([(od.accel(r + h * e) - od.accel(r - h * e)) / (2 * h)
                               for e in np.eye(3)])
        assert np.allclose(G, Gfd, rtol=1e-7, atol=1e-14)
        assert np.allclose(G, G.T, atol=1e-16)       # Hessian of a potential


def test_measurement_partials_vs_finite_difference():
    R, V = od.station_eci(od.geodetic_to_ecef(13.0, 77.5, 0.9), 1234.0)
    x = od.X_TRUE.copy()
    Ht = od.meas_partials(x, R, V)
    Hfd = np.empty((2, 6))
    for k in range(6):
        h = 1e-4 if k < 3 else 1e-7
        e = np.zeros(6); e[k] = h
        Hfd[:, k] = (np.array(od.range_rangerate(x + e, R, V))
                     - np.array(od.range_rangerate(x - e, R, V))) / (2 * h)
    assert np.allclose(Ht, Hfd, rtol=1e-6, atol=1e-12)


def test_stm_vs_perturbed_propagation():
    t = np.array([0.0, 3 * 3600.0])
    X, Phi = od.propagate(od.X_TRUE, t)
    for k in range(6):
        h = 1e-3 if k < 3 else 1e-6               # 1 m, 1 mm/s
        e = np.zeros(6); e[k] = h
        Xp, _ = od.propagate(od.X_TRUE + e, t)
        Xm, _ = od.propagate(od.X_TRUE - e, t)
        col = (Xp[-1] - Xm[-1]) / (2 * h)
        assert np.allclose(Phi[-1][:, k], col, rtol=1e-5, atol=1e-6 * np.abs(col).max())
    # Hamiltonian flow: det(Phi) = 1.
    assert abs(np.linalg.det(Phi[-1]) - 1) < 1e-8


def test_noise_free_recovers_truth():
    Y = od.simulate(od.X_TRUE, OBS)
    x0 = od.X_TRUE + np.array([2.0, -1.5, 1.0, 2e-3, -1e-3, 1.5e-3])
    x, _, info = od.batch_ls(x0, OBS, Y)
    err = x - od.X_TRUE
    assert np.linalg.norm(err[:3]) < 1e-6           # < 1 mm
    assert np.linalg.norm(err[3:]) < 1e-9           # < 1 um/s
    assert info["cost"] < 1e-4


def test_monte_carlo_covariance_consistency():
    """Formal covariance must match the scatter: NEES ~ chi2(6) per trial."""
    trials = 60
    rng = np.random.default_rng(2026)
    Y0 = od.simulate(od.X_TRUE, OBS)
    errs, nees, red_chi2 = [], [], []
    P_ref = None
    for _ in range(trials):
        Y = Y0 + rng.standard_normal(Y0.shape) * [od.SIG_RHO, od.SIG_RDOT]
        x0 = od.X_TRUE + np.concatenate([rng.normal(0, 0.1, 3), rng.normal(0, 1e-4, 3)])
        x, P, info = od.batch_ls(x0, OBS, Y)
        e = x - od.X_TRUE
        errs.append(e)
        nees.append(e @ np.linalg.solve(P, e))
        red_chi2.append(info["cost"] / info["dof"])
        P_ref = P
    total = sum(nees)
    lo, hi = chi2.ppf([0.005, 0.995], 6 * trials)
    assert lo < total < hi, (total, lo, hi)
    # Per-axis sample variance vs formal variance (chi2 with `trials` dof).
    var_ratio = np.var(np.array(errs), axis=0, ddof=0) / np.diag(P_ref)
    lo1, hi1 = chi2.ppf([0.0005, 0.9995], trials) / trials
    assert np.all((var_ratio > lo1) & (var_ratio < hi1)), var_ratio
    # Post-fit cost / dof is chi2(dof)/dof, here dof ~ 890 -> 1 +/- 0.05.
    assert abs(np.mean(red_chi2) - 1) < 0.02
    print(f"  MC: mean NEES {total/trials:.2f} (6), var ratios {np.round(var_ratio, 2)}, "
          f"mean reduced chi2 {np.mean(red_chi2):.4f}")


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"PASS {t.__name__}")
    print(f"{len(tests)} tests passed")
