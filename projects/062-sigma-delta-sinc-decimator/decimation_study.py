# Astron Built this AP
# Project 062 | Signal processing / Embedded firmware | Sinc^N decimation of a second-order sigma-delta stream
# Author: Aditya (aditya-astron)
# Date: 2026-11-24 IST
"""
Sinc^N (CIC) decimation of an isolated second-order sigma-delta current-sense
stream, as used on a PCDU shunt channel.

The modulator and the fixed-point decimator are the C firmware in cic.c,
loaded here through ctypes. This script answers three design questions:

  1. How does output SNR scale with the decimation ratio R for sinc order
     N = 1..4, and does a linear noise model predict it?
  2. What is the smallest R (highest output rate, lowest latency) that meets
     80 dB SNR for a half-scale sine, for each N?
  3. At the chosen design, how many register bits can Hogenauer pruning drop,
     and what does it cost in SNR?

Linear model. With e = y - s2 the modulator obeys Y = z^-2 X + (1-z^-1)^2 E
exactly, so the decimator output is CIC(x) plus CIC(NTF e). If e is white
with variance sigma_e^2, the output noise variance is
    sigma_out^2 = sigma_e^2 * sum_k g[k]^2 / (RM)^(2N),
    G(z) = (1 - z^-1)^2 * ((1 - z^-RM) / (1 - z^-1))^N,
and decimation does not change it. For M = 1:
    N = 1: G = (1 - z^-1)(1 - z^-R)          sum g^2 = 4
    N = 2: G = (1 - z^-R)^2                  sum g^2 = 6
    N = 3: G = (1 - z^-R)^2 (1 + ... z^-(R-1)) sum g^2 = 6R
so SNR grows as R^2, R^4, R^5 (6, 12, 15 dB per octave): with N <= L the
shaped noise near fs/2 folds straight into the output band. N = L + 1 = 3
is the first order that captures the full (2L + 1) 10log10(2) per octave.
"""
import ctypes
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FS = 10e6          # modulator clock [Hz]
AMP = 0.5          # half-scale test sine (modulator full scale = +/-1)
F0 = 1.0e3         # test tone [Hz]
SNR_REQ_DB = 80.0
N_OUT = 8192       # decimated samples per SNR measurement


def _load():
    path = os.path.join(HERE, "libcic.so")
    if not os.path.exists(path):
        raise SystemExit("libcic.so not found: run `make` first")
    lib = ctypes.CDLL(path)
    i8p = np.ctypeslib.ndpointer(np.int8, flags="C")
    f8p = np.ctypeslib.ndpointer(np.float64, flags="C")
    i64p = np.ctypeslib.ndpointer(np.int64, flags="C")
    i4p = np.ctypeslib.ndpointer(np.int32, flags="C")
    lib.sdm2_run.argtypes = [f8p, ctypes.c_size_t, i8p, ctypes.c_void_p]
    lib.sdm2_run.restype = None
    lib.cic_full_width.argtypes = [ctypes.c_int] * 4
    lib.cic_full_width.restype = ctypes.c_int
    lib.cic_init.argtypes = [ctypes.c_void_p] + [ctypes.c_int] * 4 + [ctypes.c_void_p]
    lib.cic_init.restype = ctypes.c_int
    lib.cic_run_i8.argtypes = [ctypes.c_void_p, i8p, ctypes.c_size_t, i64p]
    lib.cic_run_i8.restype = ctypes.c_size_t
    lib.cic_run.argtypes = [ctypes.c_void_p, i4p, ctypes.c_size_t, i64p]
    lib.cic_run.restype = ctypes.c_size_t
    lib.cic_hogenauer_F2.argtypes = [ctypes.c_int] * 4
    lib.cic_hogenauer_F2.restype = ctypes.c_double
    lib.cic_hogenauer_prune.argtypes = [ctypes.c_int] * 5 + [ctypes.c_void_p]
    lib.cic_hogenauer_prune.restype = None
    return lib


LIB = _load()
CIC_STATE_BYTES = 1024   # generous upper bound on sizeof(cic_t)


def modulate(x):
    x = np.ascontiguousarray(x, dtype=np.float64)
    y = np.empty(x.size, dtype=np.int8)
    e = np.empty(x.size, dtype=np.float64)
    LIB.sdm2_run(x, x.size, y, e.ctypes.data_as(ctypes.c_void_p))
    return y, e


def full_width(N, R, M=1, bin_bits=2):
    return LIB.cic_full_width(N, R, M, bin_bits)


def hogenauer_prune(N, R, bout, M=1, bin_bits=2):
    w = full_width(N, R, M, bin_bits)
    p = (ctypes.c_int * (2 * N + 2))()
    LIB.cic_hogenauer_prune(N, R, M, w, bout, p)
    return w, list(p)


def cic(bits, N, R, M=1, prune=None, width=None):
    """Run the C decimator on a +/-1 int8 stream (or int32 samples).
    Returns outputs scaled back to full-precision units."""
    width = width or full_width(N, R, M)
    state = ctypes.create_string_buffer(CIC_STATE_BYTES)
    pp = None
    if prune is not None:
        pp = (ctypes.c_int * len(prune))(*prune)
    if LIB.cic_init(state, N, R, M, width, pp) != 0:
        raise ValueError("bad CIC parameters")
    out = np.empty(bits.size // R + 1, dtype=np.int64)
    if bits.dtype == np.int8:
        m = LIB.cic_run_i8(state, np.ascontiguousarray(bits), bits.size, out)
    else:
        m = LIB.cic_run(state, np.ascontiguousarray(bits, dtype=np.int32), bits.size, out)
    shift = prune[2 * N + 1] if prune is not None else 0
    return out[:m].astype(np.float64) * 2.0 ** shift


def noise_gain(N, R, M=1):
    """sum g[k]^2 for G(z) = (1-z^-1)^2 * boxcar_RM^N."""
    g = np.array([1.0, -2.0, 1.0])
    box = np.ones(R * M)
    for _ in range(N):
        g = np.convolve(g, box)
    return float(np.sum(g * g))


def droop(f, N, R, M=1, fs=FS):
    """|H(f)| / (RM)^N of the sinc^N filter."""
    x = np.pi * f / fs
    return abs(np.sin(R * M * x) / (R * M * np.sin(x))) ** N


def sine_fit_snr(y, f, fs_out):
    """Least-squares fit of a sine at known frequency plus DC; SNR in dB."""
    n = np.arange(y.size)
    A = np.column_stack([np.cos(2 * np.pi * f / fs_out * n),
                         np.sin(2 * np.pi * f / fs_out * n), np.ones(y.size)])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A @ c
    ps = 0.5 * (c[0] ** 2 + c[1] ** 2)
    return 10 * np.log10(ps / np.mean(r * r)), r


def run_tone(N, R, n_out=N_OUT, prune=None, amp=AMP, f0=F0):
    """Simulate modulator + decimator; return measured and predicted SNR."""
    skip = N + 2
    n = (n_out + skip) * R
    t = np.arange(n) / FS
    bits, e = modulate(amp * np.sin(2 * np.pi * f0 * t))
    y = cic(bits, N, R, prune=prune)[skip:] / float(R) ** N
    snr, resid = sine_fit_snr(y, f0, FS / R)
    var_e = float(np.var(e))
    p_sig = 0.5 * (amp * droop(f0, N, R)) ** 2
    p_noise = var_e * noise_gain(N, R) / float(R) ** (2 * N)
    return snr, 10 * np.log10(p_sig / p_noise), var_e


def predicted_snr(N, R, var_e, amp=AMP, f0=F0):
    return 10 * np.log10(0.5 * (amp * droop(f0, N, R)) ** 2 * R ** (2 * N)
                         / (var_e * noise_gain(N, R)))


def f3db(N, R, M=1, fs=FS):
    """-3 dB frequency of the sinc^N response (bisection on the main lobe)."""
    lo, hi = 1e-9, fs / (R * M)
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if droop(mid, N, R, M, fs) > 10 ** (-3.0 / 20):
            lo = mid
        else:
            hi = mid
    return lo


def min_ratio(N, var_e, req=SNR_REQ_DB, rmax=4096):
    """Smallest integer R whose predicted SNR meets req."""
    for R in range(2, rmax + 1):
        if predicted_snr(N, R, var_e) >= req:
            return R
    return None


def min_ratio_verified(N, var_e, req=SNR_REQ_DB):
    """Start from the white-noise prediction, then step R up until the
    simulated chain meets req. The 1-bit error is not quite white, so the
    model is about 1 dB optimistic for N >= 3."""
    R0 = min_ratio(N, var_e, req)
    if R0 is None:
        return None, None
    R = R0
    while run_tone(N, R)[0] < req:
        R += 1
    return R0, R


def pruning_noise_db(N, R, prune):
    """Hogenauer output error variance (exact discrete truncation moments),
    relative to the full-scale unit (RM)^N."""
    v = 0.0
    for j in range(1, 2 * N + 2):
        a, b = 2.0 ** prune[j], 2.0 ** prune[j - 1]
        v += (a * a - b * b) / 12.0 * LIB.cic_hogenauer_F2(N, R, 1, j)
    return v / float(R) ** (2 * N)


def main():
    print("Second-order 1-bit modulator at %.0f MHz, %.1f-scale %.0f kHz sine\n"
          % (FS / 1e6, AMP, F0 / 1e3))

    # 1. SNR vs R for each sinc order
    ratios = [16, 32, 64, 128, 256]
    print("1. SNR vs decimation ratio: simulation / linear-model prediction [dB]")
    print("   " + "R".rjust(6) + "".join(("N=%d" % N).rjust(18) for N in range(1, 5)))
    table = {}
    var_e_all = []
    for R in ratios:
        row = []
        for N in range(1, 5):
            snr, pred, var_e = run_tone(N, R)
            table[N, R] = snr
            var_e_all.append(var_e)
            row.append("%7.1f / %-7.1f" % (snr, pred))
        print("   %6d  %s" % (R, " ".join(row)))
    slopes = []
    for N in range(1, 5):
        s = np.polyfit(np.log2(ratios), [table[N, R] for R in ratios], 1)[0]
        slopes.append(s)
    print("   slope dB/octave:   " + "".join("%18.1f" % s for s in slopes))
    print("   theory:            " + "".join("%18.1f" % s for s in
                                             [6.02, 12.04, 15.05, 15.05]))
    var_e = float(np.mean(var_e_all))
    print("   quantizer error variance sigma_e^2 = %.3f (uniform-error model 1/3 for +/-1 levels)\n"
          % var_e)

    # 2. smallest R that meets the requirement
    print("2. Smallest R for SNR >= %.0f dB (half-scale sine)" % SNR_REQ_DB)
    print("   %3s %8s %6s %12s %10s %8s %10s %12s %6s" % ("N", "R model", "R sim", "out rate",
                                                         "pred dB", "sim dB", "f-3dB",
                                                         "delay [us]", "Bmax"))
    design = {}
    for N in range(1, 5):
        R0, R = min_ratio_verified(N, var_e)
        if R is None:
            print("   %3d  not reachable below R = 4096" % N)
            continue
        snr, pred, _ = run_tone(N, R)
        delay = N * (R - 1) / 2 / FS * 1e6
        design[N] = R
        print("   %3d %8d %6d %7.1f kS/s %10.1f %8.1f %7.1f kHz %9.1f %6d"
              % (N, R0, R, FS / R / 1e3, pred, snr, f3db(N, R) / 1e3, delay,
                 full_width(N, R)))
    print()

    # 3. pruning at the N = 3 design
    N = 3
    R = design[N]
    w = full_width(N, R)
    base, _, _ = run_tone(N, R)
    # measured modulator noise at the output, relative to (RM)^N
    noise_db = 10 * np.log10(0.5 * (AMP * droop(F0, N, R)) ** 2) - base
    print("3. Hogenauer pruning at N = %d, R = %d (full width %d bits)" % (N, R, w))
    print("   modulator noise at output: %.1f dB re (RM)^N" % noise_db)
    print("   %5s %-28s %14s %10s %10s" % ("bout", "discarded B_1..B_7", "prune noise",
                                           "pred loss", "sim loss"))
    for bout in (w, 18, 17, 16, 15, 14, 12):
        _, prune = hogenauer_prune(N, R, bout)
        pn = pruning_noise_db(N, R, prune)
        snr, _, _ = run_tone(N, R, prune=prune)
        loss_pred = 10 * np.log10(1 + pn / 10 ** (noise_db / 10)) if pn > 0 else 0.0
        print("   %5d %-28s %11.1f dB %10.2f %10.2f"
              % (bout, " ".join(str(b) for b in prune[1:]),
                 10 * np.log10(pn) if pn > 0 else -np.inf, loss_pred, base - snr))

    _, prune = hogenauer_prune(N, R, 18)
    snr, _, _ = run_tone(N, R, prune=prune)
    print("\nDesign: sinc^%d, R = %d -> %.1f kS/s, %.1f us group delay, 18-bit output word,"
          % (N, R, FS / R / 1e3, N * (R - 1) / 2 / FS * 1e6))
    print("        registers %s bits, simulated SNR %.1f dB"
          % ("/".join(str(w - b) for b in prune[1:2 * N + 1]), snr))
    return 0


if __name__ == "__main__":
    sys.exit(main())
