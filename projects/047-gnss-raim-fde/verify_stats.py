# Astron Built this AP
# Project 047 | Navigation | GNSS snapshot positioning with RAIM fault detection and exclusion
# Author: Aditya (aditya-astron)
# Date: 2026-11-19 IST
"""Cross-check the C++ chi-square quantiles and noncentral CDFs against scipy.

Usage: python3 verify_stats.py ./test_raim
"""

import subprocess
import sys

from scipy.stats import chi2, ncx2


def main(binary: str) -> int:
    out = subprocess.run([binary, "--dump"], capture_output=True, text=True, check=True).stdout
    worst_inv = worst_ncx = 0.0
    n = 0
    for line in out.splitlines():
        kind, *vals = line.split()
        if kind == "inv":
            k, p, got = int(vals[0]), float(vals[1]), float(vals[2])
            # isf on the upper tail keeps scipy accurate as p -> 1.
            want = chi2.isf(1.0 - p, k) if p > 0.5 else chi2.ppf(p, k)
            worst_inv = max(worst_inv, abs(got / want - 1.0))
        else:
            k, lam, x, got = int(vals[0]), float(vals[1]), float(vals[2]), float(vals[3])
            want = ncx2.cdf(x, k, lam)
            worst_ncx = max(worst_ncx, abs(got - want))
        n += 1
    print(f"{n} values; worst chi2_inv relative error {worst_inv:.2e}, worst ncchi2_cdf absolute error {worst_ncx:.2e}")
    ok = worst_inv < 1e-9 and worst_ncx < 1e-10
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "./test_raim"))
