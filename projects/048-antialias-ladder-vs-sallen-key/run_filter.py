# Astron Built this AP
# Project 048 | Electronics and circuit analysis | ADC anti-alias filter: Chebyshev LC ladder vs. Sallen-Key cascade
# Author: Aditya (aditya-astron)
# Date: 2026-11-20 IST
"""
Anti-alias filter for a 40 MS/s telemetry ADC fed from a 50 ohm IF chain.

Spec: 0.5 dB ripple to 8 MHz, >= 60 dB at 32 MHz (fs - f_pass, the first
frequency that folds into the band). Production acceptance allows 1.0 dB
peak-to-peak passband ripple and 60 dB rejection after tolerances. Ripple is
max - min gain over 0..8 MHz, so an active stage that peaks above 0 dB is
penalised the same as one that sags.

1. Order: Chebyshev vs. Butterworth.
2. Synthesise the LC ladder and the Sallen-Key cascade on the same poles.
3. Ladder dissipation from finite inductor Q, predicted from phase
   sensitivities and checked against the full lossy analysis.
4. Orchard's argument: first-order passband sensitivity of the ladder
   vanishes at every reflection zero; the active cascade has no such zeros.
5. Op-amp gain-bandwidth needed by the cascade.
6. Tolerance Monte Carlo and production yield for both.
"""

import math

import numpy as np

import filters as fl

R0 = 50.0
F_PASS = 8e6
F_STOP = 32e6
AP_DB = 0.5
AS_DB = 60.0
ACCEPT_RIPPLE_DB = 1.0
R_SK = 200.0          # keeps the smallest Sallen-Key capacitor near 10 pF
Q_L = 60.0            # 0805 wirewound chip inductor around 8 MHz
Q_C = 1000.0          # C0G/NP0
GBW_FAST = 4e9        # fast voltage-feedback op-amp for the "realistic" row
EDGE_LOSS_DB = 0.25   # dissipation budget at the band edge

F_PB = np.linspace(0.0, F_PASS, 801)
F_SB = np.array([F_STOP, 36e6, 40e6, 60e6])


def db(x):
    with np.errstate(divide="ignore"):   # |S11| = 0 exactly at DC
        return 20.0 * np.log10(np.abs(x))


def nominal_designs():
    n = fl.cheby_order(AP_DB, AS_DB, F_STOP / F_PASS)
    g = fl.cheby_g(n, AP_DB)
    elems, rl = fl.ladder_values(g, R0, F_PASS)
    poles = fl.cheby_poles(n, AP_DB) * fl.TWO_PI * F_PASS
    secs = fl.sallen_key_design(poles, R_SK)
    return n, g, elems, rl, poles, secs


def ladder_phase_sens(f, kinds, vals, rl, h=1e-6):
    """d(phase S21)/d(ln x_k) by central difference, shape (k, F)."""
    out = []
    for k in range(len(vals)):
        up = np.array(vals, float); up[k] *= math.exp(h)
        dn = np.array(vals, float); dn[k] *= math.exp(-h)
        s_up = fl.ladder_sparams(f, kinds, up, R0, rl)[0]
        s_dn = fl.ladder_sparams(f, kinds, dn, R0, rl)[0]
        out.append(np.angle(s_up / s_dn) / (2 * h))
    return np.array(out)


def mag_sens_db(resp, vals, h=1e-6):
    """d|H|dB / d(ln x_k) for a response callable of a value vector."""
    out = []
    for k in range(len(vals)):
        up = np.array(vals, float); up[k] *= math.exp(h)
        dn = np.array(vals, float); dn[k] *= math.exp(-h)
        out.append((db(resp(up)) - db(resp(dn))) / (2 * h))
    return np.array(out)


def sk_flatten(secs):
    """Sections <-> flat value vector, so tolerances can act per component."""
    vals = []
    for s in secs:
        vals.extend(s[1:])
    return np.array(vals, float)


def sk_rebuild(secs, vals):
    out, i = [], 0
    for s in secs:
        m = len(s) - 1
        out.append((s[0],) + tuple(vals[..., j] for j in range(i, i + m)))
        i += m
    return out


def sk_component_kinds(secs):
    kinds = []
    for s in secs:
        kinds.extend(["R", "R", "C", "C"] if s[0] == "SK" else ["R", "C"])
    return kinds


def min_gbw(secs, extra_db=0.1):
    """Smallest GBW keeping the passband within extra_db of the ideal cascade."""
    ideal = db(fl.sallen_key_response(F_PB, secs))

    def err(gbw):
        return np.max(np.abs(db(fl.sallen_key_response(F_PB, secs, gbw)) - ideal))

    lo, hi = 1e7, 1e11
    for _ in range(80):
        mid = math.sqrt(lo * hi)
        lo, hi = (mid, hi) if err(mid) > extra_db else (lo, mid)
    return hi


def required_q_l(kinds, vals, rl):
    """Inductor Q for EDGE_LOSS_DB at F_PASS, from the phase-sensitivity rule."""
    sens = ladder_phase_sens(np.array([F_PASS]), kinds, vals, rl)[:, 0]
    is_l = np.array([k == "L" for k in kinds])
    k_db = 20.0 / math.log(10.0)
    return -k_db * sens[is_l].sum() / (EDGE_LOSS_DB + k_db * sens[~is_l].sum() / Q_C)


def monte_carlo(n_trials, seed=48):
    n, g, elems, rl, poles, secs = nominal_designs()
    kinds = [k for k, _ in elems]
    lad0 = np.array([v for _, v in elems])
    sk0 = sk_flatten(secs)
    sk_kinds = np.array(sk_component_kinds(secs))
    rng = np.random.default_rng(seed)

    def lad_row(tol_l, tol_c, q_l=math.inf, q_c=math.inf):
        tol = np.array([tol_l if k == "L" else tol_c for k in kinds])
        x = lad0 * (1 + tol * rng.uniform(-1, 1, (n_trials, len(lad0))))
        pb = db(fl.ladder_sparams(F_PB, kinds, x, R0, rl, q_l, q_c)[0])
        sb = db(fl.ladder_sparams(F_SB, kinds, x, R0, rl, q_l, q_c)[0])
        return pb.max(axis=1) - pb.min(axis=1), -sb.max(axis=1)

    def sk_row(tol_r, tol_c, gbw=math.inf):
        tol = np.where(sk_kinds == "R", tol_r, tol_c)
        x = sk0 * (1 + tol * rng.uniform(-1, 1, (n_trials, len(sk0))))
        s = sk_rebuild(secs, x)
        pb = db(fl.sallen_key_response(F_PB, s, gbw))
        sb = db(fl.sallen_key_response(F_SB, s, gbw))
        return pb.max(axis=1) - pb.min(axis=1), -sb.max(axis=1)

    q_need = math.ceil(required_q_l(kinds, lad0, rl))
    rows = [
        ("LC ladder, L/C +-5 %", lad_row(0.05, 0.05)),
        ("LC ladder, L +-5 %, C +-2 %", lad_row(0.05, 0.02)),
        ("LC ladder, L +-2 %, C +-2 %", lad_row(0.02, 0.02)),
        ("LC ladder, L +-5 %, C +-2 %, QL 60", lad_row(0.05, 0.02, Q_L, Q_C)),
        (f"LC ladder, L/C +-2 %, QL {q_need}", lad_row(0.02, 0.02, q_need, Q_C)),
        ("Sallen-Key, R +-1 %, C +-5 %", sk_row(0.01, 0.05)),
        ("Sallen-Key, R +-1 %, C +-2 %", sk_row(0.01, 0.02)),
        ("Sallen-Key, R +-1 %, C +-1 %", sk_row(0.01, 0.01)),
        ("Sallen-Key, R/C +-1 %, 4 GHz GBW", sk_row(0.01, 0.01, GBW_FAST)),
    ]
    return rows


def main():
    n, g, elems, rl, poles, secs = nominal_designs()
    kinds = [k for k, _ in elems]
    vals = np.array([v for _, v in elems])

    print("== Order for 0.5 dB to 8 MHz, 60 dB at 32 MHz ==")
    print(f"  Chebyshev I : n = {n}")
    print(f"  Butterworth : n = {fl.butter_order(AP_DB, AS_DB, F_STOP / F_PASS)}")

    print("\n== LC ladder (50 ohm, shunt-C first) ==")
    print("  g = " + ", ".join(f"{x:.4f}" for x in g[1:-1]) + f"   g{n+1} = {g[-1]:.4f}")
    for i, (k, v) in enumerate(elems, 1):
        unit = (v * 1e12, "pF") if k == "C" else (v * 1e6, "uH")
        print(f"  {k}{i}: {unit[0]:8.2f} {unit[1]}")
    s21, s11 = fl.ladder_sparams(F_PB, kinds, vals, R0, rl)
    s21_sb = fl.ladder_sparams(F_SB, kinds, vals, R0, rl)[0]
    print(f"  passband ripple      {-db(s21).min():.4f} dB")
    print(f"  worst return loss    {-db(s11).max():.2f} dB "
          f"(theory {-10*math.log10(1 - 10**(-AP_DB/10)):.2f})")
    print(f"  attenuation 32 MHz   {-db(s21_sb[0]):.2f} dB")

    print("\n== Sallen-Key cascade (R = 200 ohm, unity gain) ==")
    for s in secs:
        if s[0] == "SK":
            w0 = 1 / (s[1] * math.sqrt(s[3] * s[4]))
            q = 0.5 * math.sqrt(s[3] / s[4])
            print(f"  SK  f0 {w0/fl.TWO_PI/1e6:6.3f} MHz  Q {q:5.3f}  "
                  f"C1 {s[3]*1e12:7.2f} pF  C2 {s[4]*1e12:6.2f} pF  (C1/C2 = {s[3]/s[4]:.1f})")
        else:
            print(f"  RC  fc {1/(fl.TWO_PI*s[1]*s[2])/1e6:6.3f} MHz  C {s[2]*1e12:7.2f} pF")
    hsk = fl.sallen_key_response(F_PB, secs)
    print(f"  max |H_SK - H_ladder| in passband {np.max(np.abs(hsk - s21)):.2e}")

    print("\n== Ladder dissipation, QL = 60, QC = 1000 ==")
    f_chk = np.array([1e6, 4e6, 6e6, 7.6e6, 8e6])
    lossy = fl.ladder_sparams(f_chk, kinds, vals, R0, rl, Q_L, Q_C)[0]
    clean = fl.ladder_sparams(f_chk, kinds, vals, R0, rl)[0]
    sens = ladder_phase_sens(f_chk, kinds, vals, rl)
    qs = np.array([Q_L if k == "L" else Q_C for k in kinds])[:, None]
    pred = -8.685889638 * np.sum(sens / qs, axis=0)
    tau = fl.pole_group_delay(fl.TWO_PI * f_chk, poles)
    print("  f MHz   delay ns  exact dB   predicted dB")
    for i, f in enumerate(f_chk):
        print(f"  {f/1e6:5.2f}  {tau[i]*1e9:8.2f}  {db(clean[i]) - db(lossy[i]):8.4f}  {pred[i]:8.4f}")
    q_need = required_q_l(kinds, vals, rl)
    q_up = math.ceil(q_need)
    edge = (db(fl.ladder_sparams([F_PASS], kinds, vals, R0, rl)[0])
            - db(fl.ladder_sparams([F_PASS], kinds, vals, R0, rl, q_up, Q_C)[0]))[0]
    print(f"  QL for {EDGE_LOSS_DB} dB extra at 8 MHz: {q_need:.1f} -> QL {q_up} "
          f"gives {edge:.4f} dB exact")

    print("\n== First-order passband sensitivity, |dB| per 1 % component change ==")
    f_rz = np.cos((2 * np.arange(1, (n + 1) // 2 + 1) - 1) * math.pi / (2 * n)) * F_PASS
    lad_resp = lambda x, f: fl.ladder_sparams(f, kinds, x, R0, rl)[0]
    sk0 = sk_flatten(secs)
    sk_resp = lambda x, f: fl.sallen_key_response(f, sk_rebuild(secs, x))
    for f in list(f_rz) + [0.4 * F_PASS, F_PASS]:
        fa = np.array([f])
        sl = np.abs(mag_sens_db(lambda x: lad_resp(x, fa), vals)).max() * 0.01
        ss = np.abs(mag_sens_db(lambda x: sk_resp(x, fa), sk0)).max() * 0.01
        tag = " (reflection zero)" if f in f_rz else ""
        print(f"  {f/1e6:5.3f} MHz  ladder {sl:9.2e}   Sallen-Key {ss:9.2e}{tag}")

    print("\n== Op-amp for the cascade ==")
    ideal = db(fl.sallen_key_response(F_PB, secs))
    hi_q = max((s for s in secs if s[0] == "SK"), key=lambda s: s[3] / s[4])
    w0 = 1 / (hi_q[1] * math.sqrt(hi_q[3] * hi_q[4]))
    q = 0.5 * math.sqrt(hi_q[3] / hi_q[4])
    for gb in (0.5e9, 1e9, 2e9, 4e9):
        err = np.max(np.abs(db(fl.sallen_key_response(F_PB, secs, gb)) - ideal))
        print(f"  GBW {gb/1e9:3.1f} GHz: passband error {err:.3f} dB, "
              f"2 Q w0 / wt = {2*q*w0/(fl.TWO_PI*gb)*100:.2f} %")
    gbw = min_gbw(secs)
    print(f"  GBW for <= 0.1 dB passband error: {gbw/1e6:.0f} MHz")

    print("\n== Monte Carlo, 20000 trials, uniform tolerances ==")
    print(f"  {'case':38s} {'ripple p50':>10s} {'p99':>6s} {'rej min':>8s} {'yield':>7s}")
    for name, (ripple, rej) in monte_carlo(20000):
        ok = (ripple <= ACCEPT_RIPPLE_DB) & (rej >= AS_DB)
        print(f"  {name:38s} {np.median(ripple):10.3f} {np.percentile(ripple, 99):6.3f} "
              f"{rej.min():8.2f} {100*ok.mean():6.2f}%")


if __name__ == "__main__":
    main()
