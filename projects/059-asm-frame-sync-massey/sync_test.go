// Astron Built this AP
// Project 059 | Telemetry framing | Optimal ASM frame synchronization: Massey rule vs correlation
// Author: Aditya (aditya-astron)
// Date: 2026-11-23 IST

package main

import (
	"math"
	"math/bits"
	"math/rand/v2"
	"strings"
	"testing"
)

const asmString = "00011010110011111111110000011101" // 0x1ACFFC1D

func TestASMBits(t *testing.T) {
	var b strings.Builder
	for i := 0; i < ASMLen; i++ {
		b.WriteByte(byte('0' + ASMBit(i)))
	}
	if b.String() != asmString {
		t.Fatalf("ASM bits %s, want %s", b.String(), asmString)
	}
}

// Overlap mismatches by direct string comparison, independent of the
// shift-and-mask bit trick.
func TestOverlapMismatchBrute(t *testing.T) {
	for k := 1; k < ASMLen; k++ {
		m := 0
		for i := 0; i+k < ASMLen; i++ {
			if asmString[i+k] != asmString[i] {
				m++
			}
		}
		if got := OverlapMismatch(k); got != m {
			t.Errorf("k=%d: m=%d, brute %d", k, got, m)
		}
	}
	// k = 1: number of bit transitions in the marker.
	if OverlapMismatch(1) != 11 {
		t.Errorf("m_1 = %d, want 11 transitions", OverlapMismatch(1))
	}
}

func TestBER(t *testing.T) {
	// Q(sqrt 2) at Es/N0 = 0 dB.
	if got := BER(1); math.Abs(got-0.0786496035) > 1e-9 {
		t.Fatalf("BER(0 dB) = %.10f", got)
	}
}

func TestDistancePMFMoments(t *testing.T) {
	// Random window, E = 2: (1 + 32 + 496)/2^32, the figure used in project 002.
	pr := DistancePMF(ASMLen, 0.2)
	pf := pr[0] + pr[1] + pr[2]
	if math.Abs(pf-529/math.Pow(2, 32)) > 1e-20 {
		t.Fatalf("random P_F(E=2) = %.6e", pf)
	}
	for _, p := range []float64{0, 0.05, 0.3} {
		for k := 0; k <= ASMLen; k++ {
			pmf := DistancePMF(k, p)
			sum, mean := 0.0, 0.0
			for d, v := range pmf {
				sum += v
				mean += float64(d) * v
			}
			kk := min(k, ASMLen)
			m := OverlapMismatch(k)
			want := float64(ASMLen-kk-m)*p + float64(m)*(1-p) + float64(kk)/2
			if math.Abs(sum-1) > 1e-12 || math.Abs(mean-want) > 1e-10 {
				t.Errorf("k=%d p=%.2f: sum %.15f mean %.6f want %.6f", k, p, sum, mean, want)
			}
		}
	}
}

// Shifted windows on a real noisy stream vs the exact convolution pmf.
func TestDistancePMFMonteCarlo(t *testing.T) {
	rng := rand.New(rand.NewPCG(59, 1))
	const n = 200000
	const p = 0.1
	for _, k := range []int{1, 7, 24} {
		counts := make([]int, ASMLen+1)
		for trial := 0; trial < n; trial++ {
			// bits: ASM followed by k random data bits, then channel flips
			var stream uint64 = uint64(ASM) << 32
			stream |= (rng.Uint64() & ((1 << uint(k)) - 1)) << uint(32-k)
			for i := 0; i < 64; i++ {
				if rng.Float64() < p {
					stream ^= 1 << uint(i)
				}
			}
			win := uint32(stream >> uint(32-k))
			counts[bits.OnesCount32(win^ASM)]++
		}
		pmf := DistancePMF(k, p)
		for d, c := range counts {
			q := 0.0
			if d < len(pmf) {
				q = pmf[d]
			}
			sd := math.Sqrt(n * q * (1 - q))
			if math.Abs(float64(c)-n*q) > 5*sd+3 {
				t.Errorf("k=%d d=%d: count %d, expected %.1f +- %.1f", k, d, c, n*q, sd)
			}
		}
	}
}

func TestMasseyIsLLR(t *testing.T) {
	rng := rand.New(rand.NewPCG(59, 2))
	for _, esn0 := range []float64{0.2, 0.5, 2} {
		lam := 2 * esn0
		for trial := 0; trial < 50; trial++ {
			y := make([]float64, ASMLen)
			c, pen := 0.0, 0.0
			for i := range y {
				y[i] = float64(1-2*int(rng.Uint64()&1)) + rng.NormFloat64()/math.Sqrt(2*esn0)
				c += ASMSym(i) * y[i]
				pen += LnCosh(lam*y[i]) / lam
			}
			got := lam * (c - pen)
			want := MasseyLLR(y, esn0)
			if math.Abs(got-want) > 1e-9*math.Max(1, math.Abs(want)) {
				t.Fatalf("Es/N0=%.1f: lam*M = %.12f, LLR = %.12f", esn0, got, want)
			}
		}
	}
}

func TestHighSNRLimit(t *testing.T) {
	rng := rand.New(rand.NewPCG(59, 3))
	y := make([]float64, ASMLen)
	c, sabs, mis := 0.0, 0.0, 0.0
	for i := range y {
		y[i] = 2 * rng.NormFloat64()
		c += ASMSym(i) * y[i]
		sabs += math.Abs(y[i])
		if (y[i] >= 0) != (ASMSym(i) > 0) {
			mis += math.Abs(y[i])
		}
	}
	if math.Abs((c-sabs)-(-2*mis)) > 1e-12 {
		t.Fatalf("c - sum|y| = %.12f, -2 sum_mis |y| = %.12f", c-sabs, -2*mis)
	}
	// lncosh(lam y)/lam = |y| - ln2/lam + O(e^{-2 lam |y|}/lam): the
	// penalty differs from sum|y| by at most 32 ln2 / lam.
	for _, lam := range []float64{1, 10, 1000} {
		pen := 0.0
		for _, v := range y {
			pen += LnCosh(lam*v) / lam
		}
		if d := sabs - pen; d < -1e-12 || d > 32*math.Ln2/lam+1e-12 {
			t.Errorf("lam=%g: sum|y| - penalty = %.3e, bound %.3e", lam, d, 32*math.Ln2/lam)
		}
	}
}

// Full CADU acquisition: Monte Carlo vs the exact-overlap independence
// approximation for hard and corr, and the Massey rule ahead of both.
func TestAcquisitionAgainstAnalytic(t *testing.T) {
	if testing.Short() {
		t.Skip("Monte Carlo acquisition")
	}
	const L = caduBits
	const trials = 800
	for _, db := range []float64{-4, 0} {
		esn0 := lin(db)
		fails := monteCarlo(L, esn0, trials, 0xACE+uint64(db+10))
		check := func(name string, mcFails int, analytic float64) {
			p := float64(mcFails) / trials
			sd := math.Sqrt(analytic * (1 - analytic) / trials)
			if math.Abs(p-analytic) > 4*sd+1.0/trials {
				t.Errorf("%+.0f dB %s: MC %.4f, analytic %.4f +- %.4f", db, name, p, analytic, sd)
			}
			t.Logf("%+.0f dB %-5s MC %.4f analytic %.4f", db, name, p, analytic)
		}
		check("hard", fails[RuleHard], 1-HardAcqAnalytic(L, BER(esn0)))
		check("corr", fails[RuleCorr], 1-CorrAcqAnalytic(L, esn0))
		if !(fails[RuleMassey] < fails[RuleCorr] && fails[RuleMassey] < fails[RuleHard]) {
			t.Errorf("%+.0f dB: Massey %d fails not below corr %d and hard %d",
				db, fails[RuleMassey], fails[RuleCorr], fails[RuleHard])
		}
		t.Logf("%+.0f dB massey MC %.4f", db, float64(fails[RuleMassey])/trials)
	}
}
