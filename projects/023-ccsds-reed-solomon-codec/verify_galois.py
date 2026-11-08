# Astron Built this AP
# Project 023 | Error correction | CCSDS Reed-Solomon (255,223) codec with errors-and-erasures decoding
# Author: Aditya (aditya-astron)
# Date: 2026-11-08 IST
#
# Independent cross-check of the C codec against the `galois` library
# (pip install galois), which builds the same code from its definition:
# GF(2^8) mod x^8+x^7+x^2+x+1, primitive element beta = alpha^11, first
# consecutive root beta^112. Loads librs.so through ctypes.

import ctypes
import os
import sys

import galois
import numpy as np

N, K = 255, 223

here = os.path.dirname(os.path.abspath(__file__))
lib = ctypes.CDLL(os.path.join(here, "librs.so"))
lib.rs_init()
Buf = ctypes.c_uint8 * N

GF = galois.GF(2**8, irreducible_poly=0x187, primitive_element=2)
beta = GF(2) ** 11
rs = galois.ReedSolomon(N, K, field=GF, alpha=beta, c=112)

failures = 0


def check(cond, msg):
    global failures
    if not cond:
        failures += 1
        print("FAIL:", msg)


genpoly = np.array((ctypes.c_uint8 * 33).in_dll(lib, "rs_genpoly"))
check(np.array_equal(genpoly, np.array(rs.generator_poly.coeffs)), "generator polynomial differs")

rng = np.random.default_rng(2026)
for trial in range(50):
    msg = rng.integers(0, 256, K, dtype=np.uint8)
    buf = Buf(*msg.tolist(), *([0] * (N - K)))
    lib.rs_encode_conv(buf)
    ours = np.array(buf)
    ref = np.array(rs.encode(GF(msg)))
    check(np.array_equal(ours, ref), f"encoding differs on trial {trial}")

    # corrupt with 16 errors and let both decoders repair it
    rx = ref.copy()
    pos = rng.choice(N, 16, replace=False)
    rx[pos] ^= rng.integers(1, 256, 16, dtype=np.uint8)
    buf = Buf(*rx.tolist())
    n_ours = lib.rs_decode_conv(buf, None, 0, 0)
    dec_ref, n_ref = rs.decode(GF(rx), output="codeword", errors=True)
    check(n_ours == n_ref == 16, f"error counts {n_ours} vs {n_ref}")
    check(np.array_equal(np.array(buf), np.array(dec_ref)), f"decoded word differs on trial {trial}")

print("galois cross-check:", "OK" if failures == 0 else f"{failures} failures")
sys.exit(1 if failures else 0)
