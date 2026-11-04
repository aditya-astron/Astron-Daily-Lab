# Astron Built this AP
# Project 013 | Navigation / Signal processing | GPS L1 C/A code generation and FFT acquisition with CFAR detection
# Author: Aditya (aditya-astron)
# Date: 2026-11-04 IST
"""GPS L1 C/A Gold code generator (IS-GPS-200).

Two 10-stage maximal-length LFSRs clocked at 1.023 MHz:

    G1: 1 + x^3 + x^10
    G2: 1 + x^2 + x^3 + x^6 + x^8 + x^9 + x^10

Both start all-ones. G1 output is stage 10. The G2 contribution for a given
PRN is the modulo-2 sum of two selected stages (the "phase selector"), which
is equivalent to a fixed delay of the G2 sequence. Different delays of G2
summed with G1 give the members of a Gold family, so cross-correlation
between any two codes takes only the values {-65, -1, 63} out of 1023.
"""
import numpy as np

CODE_LENGTH = 1023
CHIP_RATE = 1.023e6  # chips/s

# G2 phase-selector taps for PRN 1..32 (IS-GPS-200, Table 3-Ia).
G2_TAPS = {
    1: (2, 6), 2: (3, 7), 3: (4, 8), 4: (5, 9), 5: (1, 9), 6: (2, 10),
    7: (1, 8), 8: (2, 9), 9: (3, 10), 10: (2, 3), 11: (3, 4), 12: (5, 6),
    13: (6, 7), 14: (7, 8), 15: (8, 9), 16: (9, 10), 17: (1, 4), 18: (2, 5),
    19: (3, 6), 20: (4, 7), 21: (5, 8), 22: (6, 9), 23: (1, 3), 24: (4, 6),
    25: (5, 7), 26: (6, 8), 27: (7, 9), 28: (8, 10), 29: (1, 6), 30: (2, 7),
    31: (3, 8), 32: (4, 9),
}


def ca_code_bits(prn):
    """One period of the C/A code for `prn` as logic bits (0/1)."""
    if prn not in G2_TAPS:
        raise ValueError(f"PRN must be 1..32, got {prn}")
    t1, t2 = G2_TAPS[prn]
    g1 = [1] * 10  # g1[0] is stage 1
    g2 = [1] * 10
    out = np.empty(CODE_LENGTH, dtype=np.uint8)
    for i in range(CODE_LENGTH):
        out[i] = g1[9] ^ g2[t1 - 1] ^ g2[t2 - 1]
        fb1 = g1[2] ^ g1[9]
        fb2 = g2[1] ^ g2[2] ^ g2[5] ^ g2[7] ^ g2[8] ^ g2[9]
        g1 = [fb1] + g1[:9]
        g2 = [fb2] + g2[:9]
    return out


def ca_code(prn):
    """C/A code as bipolar chips: logic 0 -> +1, logic 1 -> -1."""
    return 1 - 2 * ca_code_bits(prn).astype(np.int8)


def sample_code(prn, fs, n_samples, code_phase_chips=0.0):
    """Sample the code at rate fs starting `code_phase_chips` into the period.

    Nearest-chip (zero-order hold) sampling, as an ADC would see a
    rectangular-chip signal through a wide front end.
    """
    chips = ca_code(prn)
    idx = np.floor(code_phase_chips + np.arange(n_samples) * CHIP_RATE / fs)
    return chips[idx.astype(np.int64) % CODE_LENGTH].astype(np.float64)
