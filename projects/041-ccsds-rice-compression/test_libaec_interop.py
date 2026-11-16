#!/usr/bin/env python3
# Astron Built this AP
# Project 041 | Telemetry data compression | CCSDS 121.0 adaptive Rice coder
# Author: Aditya (aditya-astron)
# Date: 2026-11-16 IST
"""Cross-check the C coder against libaec (the reference CCSDS 121.0
implementation, via imagecodecs): every stream we write must decode in
libaec, every libaec stream must decode in ours, and our exhaustive option
search must never be longer than libaec's heuristic search."""
import os
import subprocess
import sys
import tempfile

import numpy as np

try:
    import imagecodecs
except ImportError:
    print("SKIP: pip install imagecodecs (provides libaec)")
    sys.exit(0)

HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.join(HERE, "rice_tool")
FLAGS = imagecodecs.AEC.FLAG.DATA_PREPROCESS


def cases(rng):
    """Signals chosen so that every coding option and run rule is hit."""
    N = 64 * 16 * 40
    t = np.arange(N)
    for n in (8, 12, 14, 16):
        top = (1 << n) - 1
        mid = top // 2
        yield n, "uniform", rng.integers(0, top + 1, N)
        for sigma in (0.3, 2.0, 40.0):
            s = mid + 0.3 * mid * np.sin(2 * np.pi * t / 5000) + rng.normal(0, sigma, N)
            yield n, f"sine+noise{sigma}", np.clip(np.round(s), 0, top)
        # plateaus of every length: zero-block runs, ROS, runs at RSI ends
        lev = np.repeat(rng.integers(0, top + 1, N // 37 + 1), 37)[:N]
        lev[rng.random(N) < 0.002] ^= 1
        yield n, "plateaus", lev
        yield n, "constant", np.full(N, mid // 3)
        # rail-hugging samples exercise the one-sided branch of the mapper
        yield n, "rails", np.where(rng.random(N) < 0.5, rng.integers(0, 4, N),
                                    top - rng.integers(0, 4, N))


def main():
    rng = np.random.default_rng(121)
    fails = 0
    total = 0
    worse = 0
    identical = 0
    with tempfile.TemporaryDirectory() as tmp:
        raw = os.path.join(tmp, "x.u16")
        ours = os.path.join(tmp, "x.rice")
        back = os.path.join(tmp, "y.u16")
        for n, name, x in cases(rng):
            for J, rsi in ((16, 128), (32, 64), (64, 8), (8, 200)):
                x16 = np.asarray(x).astype(np.uint16)
                # a partial final RSI must also round-trip
                x16 = x16[: len(x16) - J * 3]
                x16.tofile(raw)
                N = len(x16)
                subprocess.run([TOOL, "enc", str(n), str(J), str(rsi), raw, ours], check=True)
                with open(ours, "rb") as f:
                    mine = f.read()
                # libaec stores n <= 8 samples in one byte, 9..16 in two
                dt = np.uint8 if n <= 8 else np.uint16
                width = np.dtype(dt).itemsize
                ref = imagecodecs.aec_encode(x16.astype(dt).tobytes(), bitspersample=n,
                                             flags=FLAGS, blocksize=J, rsi=rsi)
                # ours -> libaec
                dec1 = np.frombuffer(imagecodecs.aec_decode(
                    mine, bitspersample=n, flags=FLAGS, blocksize=J, rsi=rsi,
                    out=width * N), dtype=dt).astype(np.uint16)
                # libaec -> ours
                with open(ours, "wb") as f:
                    f.write(ref)
                subprocess.run([TOOL, "dec", str(n), str(J), str(rsi), str(N), ours, back],
                               check=True)
                dec2 = np.fromfile(back, dtype=np.uint16)
                ok = np.array_equal(dec1, x16) and np.array_equal(dec2, x16)
                total += 1
                identical += mine == ref
                if len(mine) > len(ref):
                    worse += 1
                if not ok:
                    fails += 1
                    print(f"FAIL n={n} {name} J={J} rsi={rsi}")
                if J == 16 and rsi == 128:
                    print(f"  n={n:2d} {name:16s} ours {len(mine):7d} B   libaec {len(ref):7d} B"
                          f"   raw {N * n // 8:7d} B")
    print(f"{total - fails}/{total} interop cases round-trip both ways")
    print(f"ours longer than libaec in {worse} cases, byte-identical in {identical}")
    ok = fails == 0 and worse == 0
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
