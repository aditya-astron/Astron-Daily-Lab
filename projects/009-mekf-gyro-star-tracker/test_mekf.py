# Astron Built this AP
# Project 009 | Kalman filtering and sensor fusion | Multiplicative EKF fusing a drifting gyro with a star tracker
# Author: Aditya (aditya-astron)
# Date: 2026-11-03 IST
"""Checks against closed forms, an independent Riccati iteration and Monte Carlo.

Run: python3 -m pytest -q
"""
import numpy as np
from scipy.stats import chi2

from mekf import (ARCSEC, DEG_PER_HR, discretize, farrenkopf_qd, qerr, qexp,
                  qmul, simulate, single_axis_steady_state)

SV, SU, SN = np.sqrt(10) * 1e-7, np.sqrt(10) * 1e-10, 6.0 * ARCSEC


def test_quaternion_error_roundtrip():
    rng = np.random.default_rng(3)
    for _ in range(50):
        q = qexp(rng.normal(0, 2, 3))
        v = rng.normal(0, 0.8, 3)
        assert np.allclose(qerr(q, qmul(q, qexp(v))), v, atol=1e-12)


def test_van_loan_matches_farrenkopf_at_zero_rate():
    dt = 0.1
    Phi, Qd = discretize(np.zeros(3), dt, SV, SU)
    Phi_ref = np.block([[np.eye(3), -dt * np.eye(3)], [np.zeros((3, 3)), np.eye(3)]])
    assert np.allclose(Phi, Phi_ref, atol=1e-15)
    assert np.allclose(Qd, farrenkopf_qd(dt, SV, SU), rtol=1e-9, atol=1e-30)


def test_rotating_phi_is_body_rotation():
    # With zero bias error the attitude block of Phi is exp(-[w x] dt), i.e. the
    # transpose of the body rotation over the step.
    w = np.array([0.01, -0.02, 0.03])
    dt = 5.0
    Phi, _ = discretize(w, dt, SV, SU)
    a = np.linalg.norm(w * dt)
    k = w / np.linalg.norm(w)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    R = np.eye(3) + np.sin(a) * K + (1 - np.cos(a)) * K @ K
    assert np.allclose(Phi[:3, :3], R.T, atol=1e-12)


def test_dare_matches_brute_force_riccati_with_substeps():
    # Ten 0.1 s gyro propagations per 1 s update, iterated to convergence.
    Phi, Qd = discretize(np.zeros(3), 0.1, SV, SU)
    Phi, Qd = Phi[[0, 3]][:, [0, 3]], Qd[[0, 3]][:, [0, 3]]
    H = np.array([[1.0, 0.0]])
    P = np.diag([1e-6, 1e-12])
    for _ in range(20000):
        for _ in range(10):
            P = Phi @ P @ Phi.T + Qd
        K = P @ H.T / (H @ P @ H.T + SN**2)
        P = (np.eye(2) - K @ H) @ P
    pa, pb = single_axis_steady_state(1.0, SV, SU, SN)
    assert np.isclose(P[0, 0], pa, rtol=1e-6)
    assert np.isclose(P[1, 1], pb, rtol=1e-6)


def test_no_bias_drift_limit_matches_scalar_random_walk():
    # As sigma_u -> 0 the bias becomes perfectly known and the attitude channel
    # is a random walk (q = sigma_v^2 dt) seen through noise R:
    #   P- = q/2 + sqrt(q^2/4 + q R),  P+ = P- R / (P- + R)
    q, R = SV**2 * 1.0, SN**2
    Pm = q / 2 + np.sqrt(q**2 / 4 + q * R)
    Pp = Pm * R / (Pm + R)
    pa, _ = single_axis_steady_state(1.0, SV, 1e-16, SN)
    assert np.isclose(pa, Pp, rtol=1e-3)
    # Bias drift must make things worse, not better.
    pa_drift, _ = single_axis_steady_state(1.0, SV, SU, SN)
    assert pa_drift > Pp


def test_monte_carlo_consistency_and_accuracy():
    runs, T = 12, 1800.0
    pa, pb = single_axis_steady_state(1.0, SV, SU, SN)
    nees_end, att_sq, z_att, z_bias, sig_end = [], [], [], [], []
    for s in range(runs):
        r = simulate(T=T, seed=100 + s)
        ss = r["t"] > 600.0
        nees_end.append(r["nees"][-1])
        att_sq.append(np.mean(r["err"][ss, :3] ** 2))
        # The bias channel is still settling at 600 s (sigma ~1.5x steady
        # state), so compare errors against the filter's own sigma.
        z = r["err"][ss] / r["sig"][ss]
        z_att.append(np.mean(z[:, :3] ** 2))
        z_bias.append(np.mean(z[:, 3:] ** 2))
        sig_end.append(r["sig"][-1])
        # Bias was 0.1-0.15 deg/hr; it must be learnt to within 3 sigma.
        assert np.all(np.abs(r["err"][-1, 3:]) < 3 * r["sig"][-1, 3:])

    # Ensemble NEES at the final epoch: sum ~ chi2 with 6*runs dof.
    lo, hi = chi2.ppf([0.005, 0.995], 6 * runs) / runs
    assert lo < np.mean(nees_end) < hi, (np.mean(nees_end), lo, hi)

    # Filter covariance converges to the analytic single-axis solution; the
    # slew couples axes only weakly at these rates.
    sig_end = np.array(sig_end)
    assert np.allclose(sig_end[:, :3], np.sqrt(pa), rtol=0.02)
    assert np.allclose(sig_end[:, 3:], np.sqrt(pb), rtol=0.02)

    # Actual errors match the predicted ones (normalised RMS ~ 1).
    assert np.isclose(np.sqrt(np.mean(z_att)), 1.0, atol=0.15)
    assert np.isclose(np.sqrt(np.mean(z_bias)), 1.0, atol=0.25)
    assert np.isclose(np.sqrt(np.mean(att_sq)), np.sqrt(pa), rtol=0.15)
    # And the fused estimate beats the raw star tracker by a wide margin.
    assert np.sqrt(np.mean(att_sq)) < SN / 5


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-q"]))
