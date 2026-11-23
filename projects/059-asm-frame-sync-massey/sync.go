// Astron Built this AP
// Project 059 | Telemetry framing | Optimal ASM frame synchronization: Massey rule vs correlation
// Author: Aditya (aditya-astron)
// Date: 2026-11-23 IST

// Frame synchronization on the CCSDS 32-bit attached sync marker (ASM).
//
// Channel model: BPSK symbols x = +-1 with unit energy, y = x + n,
// n ~ N(0, N0/2), so Es/N0 = 1/N0 and the hard-decision BER is
// p = Q(sqrt(2 Es/N0)). Frame data is randomized, so data symbols are
// modelled as independent equiprobable +-1.
//
// A window of 32 samples starting at candidate offset j is scored by
//
//   hard:    -d_H(sign(y), ASM)                (count of agreeing bits)
//   corr:    sum s_i y_i                        (soft correlation)
//   massey:  sum s_i y_i - (1/lam) sum lncosh(lam y_i),   lam = 2 Es/N0
//   hisnr:   sum s_i y_i - sum |y_i|  = -2 sum_{sign(y_i) != s_i} |y_i|
//
// The Massey metric is lam^-1 times the log-likelihood ratio of "ASM here"
// against "random data here": per sample p(y|s)/p(y|random) =
// exp(lam y s)/cosh(lam y). The cosh term penalizes windows that sit on
// large-amplitude samples, which plain correlation rewards whatever their
// sign. hisnr is its lam -> infinity limit and needs no SNR estimate.

package main

import (
	"math"
	"math/bits"
	"math/rand/v2"
)

// ASM is the CCSDS 131.0-B attached sync marker, transmitted MSB first.
const ASM uint32 = 0x1ACFFC1D

// ASMLen is the marker length in bits.
const ASMLen = 32

// ASMBit returns bit i of the marker (i = 0 is sent first).
func ASMBit(i int) int { return int(ASM>>(31-i)) & 1 }

// ASMSym returns the BPSK symbol of bit i: bit 0 -> +1, bit 1 -> -1.
func ASMSym(i int) float64 { return 1 - 2*float64(ASMBit(i)) }

// OverlapMismatch returns m_k, the number of disagreements between ASM[i+k]
// and ASM[i] for i = 0..31-k. A window shifted k bits from the true marker
// (either direction) overlaps 32-k marker bits, of which m_k disagree with
// the reference; the other k bits are random data.
func OverlapMismatch(k int) int {
	if k <= 0 || k >= ASMLen {
		return 0
	}
	mask := ^uint32(0) << uint(k)
	return bits.OnesCount32(((ASM << uint(k)) ^ ASM) & mask)
}

// Binom returns C(n, k) as a float64.
func Binom(n, k int) float64 {
	if k < 0 || k > n {
		return 0
	}
	r := 1.0
	for i := 1; i <= k; i++ {
		r = r * float64(n-k+i) / float64(i)
	}
	return r
}

// BinomPMF returns the pmf of Bin(n, p) on 0..n.
func BinomPMF(n int, p float64) []float64 {
	out := make([]float64, n+1)
	for k := 0; k <= n; k++ {
		out[k] = Binom(n, k) * math.Pow(p, float64(k)) * math.Pow(1-p, float64(n-k))
	}
	return out
}

func convolve(a, b []float64) []float64 {
	out := make([]float64, len(a)+len(b)-1)
	for i, x := range a {
		for j, y := range b {
			out[i+j] += x * y
		}
	}
	return out
}

// DistancePMF is the exact pmf of the Hamming distance between the hard-decision
// window and the ASM, for a window shifted k bits from the true marker
// (k = 0: aligned, k >= 32: pure random data), at channel BER p.
//
//	D = Bin(32-k-m_k, p)  agreeing overlap bits flipped by the channel
//	  + Bin(m_k, 1-p)     disagreeing overlap bits not flipped
//	  + Bin(k, 1/2)       random data bits
func DistancePMF(k int, p float64) []float64 {
	if k >= ASMLen {
		return BinomPMF(ASMLen, 0.5)
	}
	m := OverlapMismatch(k)
	d := convolve(BinomPMF(ASMLen-k-m, p), BinomPMF(m, 1-p))
	d = convolve(d, BinomPMF(k, 0.5))
	return d
}

// QFunc is the Gaussian tail probability.
func QFunc(x float64) float64 { return 0.5 * math.Erfc(x/math.Sqrt2) }

// BER returns the uncoded BPSK hard-decision bit error rate at Es/N0 (linear).
func BER(esn0 float64) float64 { return QFunc(math.Sqrt(2 * esn0)) }

// LnCosh is a stable ln(cosh(x)).
func LnCosh(x float64) float64 {
	a := math.Abs(x)
	return a + math.Log1p(math.Exp(-2*a)) - math.Ln2
}

// MasseyLLR evaluates the per-window log-likelihood ratio directly from the
// Gaussian densities. Used to check the closed-form metric.
func MasseyLLR(y []float64, esn0 float64) float64 {
	n0 := 1 / esn0
	g := func(v, mean float64) float64 { return math.Exp(-(v - mean) * (v - mean) / n0) }
	llr := 0.0
	for i := 0; i < ASMLen; i++ {
		num := g(y[i], ASMSym(i))
		den := 0.5*g(y[i], 1) + 0.5*g(y[i], -1)
		llr += math.Log(num / den)
	}
	return llr
}

// Rule identifies a sync metric.
type Rule int

const (
	RuleHard Rule = iota
	RuleCorr
	RuleMassey
	RuleMasseyHi  // lam -> infinity, SNR-free
	RuleMasseyLo3 // Massey with Es/N0 underestimated by 3 dB
	RuleMasseyHi3 // Massey with Es/N0 overestimated by 3 dB
	NumRules
)

var RuleNames = [NumRules]string{
	"hard", "corr", "massey", "massey-hiSNR", "massey(-3dB est)", "massey(+3dB est)",
}

// Searcher holds scratch buffers for one-shot acquisition over a frame.
type Searcher struct {
	L     int // frame length in bits including the ASM (candidate offsets 0..L-1)
	esn0  float64
	sigma float64
	sym   [ASMLen]float64
	y     []float64
	pen   [3][]float64 // prefix sums of the lncosh penalty for the three lam values
	lam   [3]float64
	pabs  []float64 // prefix sums of |y|
}

// NewSearcher prepares an acquisition search for frame length L at Es/N0.
func NewSearcher(L int, esn0 float64) *Searcher {
	s := &Searcher{L: L, esn0: esn0, sigma: math.Sqrt(0.5 / esn0)}
	for i := 0; i < ASMLen; i++ {
		s.sym[i] = ASMSym(i)
	}
	n := L + ASMLen - 1
	s.y = make([]float64, n)
	for r := 0; r < 3; r++ {
		s.pen[r] = make([]float64, n+1)
	}
	s.pabs = make([]float64, n+1)
	s.lam = [3]float64{2 * esn0, 2 * esn0 * 0.5, 2 * esn0 * 2}
	return s
}

// Trial transmits ASM + random data + the first 31 bits of the next ASM,
// scores every candidate offset 0..L-1 with every rule and reports which
// rules picked the true offset 0.
func (s *Searcher) Trial(rng *rand.Rand) (correct [NumRules]bool) {
	n := len(s.y)
	for i := 0; i < n; i++ {
		var x float64
		switch {
		case i < ASMLen:
			x = s.sym[i]
		case i >= s.L:
			x = s.sym[i-s.L]
		default:
			x = float64(1 - 2*int(rng.Uint64()&1))
		}
		s.y[i] = x + s.sigma*rng.NormFloat64()
	}
	for i := 0; i < n; i++ {
		s.pabs[i+1] = s.pabs[i] + math.Abs(s.y[i])
		for r := 0; r < 3; r++ {
			s.pen[r][i+1] = s.pen[r][i] + LnCosh(s.lam[r]*s.y[i])/s.lam[r]
		}
	}

	best := [NumRules]float64{}
	arg := [NumRules]int{}
	ties := 0
	for r := range best {
		best[r] = math.Inf(-1)
	}
	// Hard decisions as a sliding 32-bit register: bit 1 where y < 0, which
	// matches the ASM bit mapping, so d_H = popcount(reg XOR ASM).
	var reg uint32
	for i := 0; i < ASMLen-1; i++ {
		reg = reg<<1 | b2u(s.y[i] < 0)
	}
	for j := 0; j < s.L; j++ {
		reg = reg<<1 | b2u(s.y[j+ASMLen-1] < 0)
		agree := ASMLen - bits.OnesCount32(reg^ASM)
		c := 0.0
		w := s.y[j : j+ASMLen]
		for i, v := range w {
			c += s.sym[i] * v
		}
		// Hard decisions tie often; break ties uniformly (reservoir rule).
		h := float64(agree)
		if h > best[RuleHard] {
			best[RuleHard], arg[RuleHard], ties = h, j, 1
		} else if h == best[RuleHard] {
			ties++
			if rng.IntN(ties) == 0 {
				arg[RuleHard] = j
			}
		}
		score := [NumRules]float64{
			RuleCorr:      c,
			RuleMassey:    c - (s.pen[0][j+ASMLen] - s.pen[0][j]),
			RuleMasseyHi:  c - (s.pabs[j+ASMLen] - s.pabs[j]),
			RuleMasseyLo3: c - (s.pen[1][j+ASMLen] - s.pen[1][j]),
			RuleMasseyHi3: c - (s.pen[2][j+ASMLen] - s.pen[2][j]),
		}
		for r := RuleCorr; r < NumRules; r++ {
			if score[r] > best[r] {
				best[r], arg[r] = score[r], j
			}
		}
	}
	for r := range arg {
		correct[r] = arg[r] == 0
	}
	return
}

func b2u(b bool) uint32 {
	if b {
		return 1
	}
	return 0
}

// candidateShifts lists, for the L candidate offsets other than the true
// one, how far each is from a marker: offsets 1..31 overlap the tail of this
// marker, offsets L-31..L-1 overlap the head of the next, the rest see
// only data. Each overlap shift k occurs twice (once per marker), so the
// returned list holds k once and callers square its factor.
func candidateShifts(L int) (overlap []int, nRandom int) {
	for k := 1; k < ASMLen; k++ {
		overlap = append(overlap, k)
	}
	return overlap, L - 1 - 2*(ASMLen-1)
}

// HardAcqAnalytic is the probability that minimum-distance search with
// random tie-breaking picks the true offset, treating the L window distances
// as independent (they share no symbols except near the marker, where the
// exact overlap pmfs are used). With D_T = d, winning against K tied windows
// has probability E[1/(1+K)] = int_0^1 prod_j (P(D_j > d) + P(D_j = d) u) du.
func HardAcqAnalytic(L int, p float64) float64 {
	pt := DistancePMF(0, p)
	overlap, nRand := candidateShifts(L)
	pmfs := [][]float64{}
	for _, k := range overlap {
		pmfs = append(pmfs, DistancePMF(k, p))
	}
	pr := DistancePMF(ASMLen, p)
	tail := func(pmf []float64, d int) (gt, eq float64) {
		for i := d + 1; i < len(pmf); i++ {
			gt += pmf[i]
		}
		return gt, pmf[d]
	}
	total := 0.0
	const nu = 4000 // composite Simpson panels on u in [0, 1]
	for d := 0; d <= ASMLen; d++ {
		if pt[d] < 1e-300 {
			continue
		}
		gr, er := tail(pr, d)
		var g, e []float64
		for _, pmf := range pmfs {
			a, b := tail(pmf, d)
			g, e = append(g, a), append(e, b)
		}
		f := func(u float64) float64 {
			lg := float64(nRand) * math.Log(gr+er*u)
			for i := range g {
				lg += 2 * math.Log(g[i]+e[i]*u)
			}
			return math.Exp(lg)
		}
		h := 1.0 / nu
		sum := f(0) + f(1)
		for i := 1; i < nu; i++ {
			w := 2.0
			if i%2 == 1 {
				w = 4
			}
			sum += w * f(float64(i)*h)
		}
		total += pt[d] * sum * h / 3
	}
	return total
}

// CorrAcqAnalytic is the same independence approximation for soft
// correlation. A window shifted k bits has correlation
// (32 - 2 m_k - 2B) + N(0, 32 sigma^2), B ~ Bin(k, 1/2); the true window has
// mean 32. P = int f_T(t) prod_j F_j(t) dt, integrated on a fine grid.
func CorrAcqAnalytic(L int, esn0 float64) float64 {
	sd := math.Sqrt(32 * 0.5 / esn0)
	overlap, nRand := candidateShifts(L)
	type comp struct{ mean, w []float64 }
	mk := func(k int) comp {
		var c comp
		m := 0
		if k < ASMLen {
			m = OverlapMismatch(k)
		} else {
			k = ASMLen
		}
		for b, pb := range BinomPMF(k, 0.5) {
			c.mean = append(c.mean, float64(ASMLen-2*m-2*b))
			c.w = append(c.w, pb)
		}
		return c
	}
	cdf := func(c comp, t float64) float64 {
		s := 0.0
		for i := range c.mean {
			s += c.w[i] * (1 - QFunc((t-c.mean[i])/sd))
		}
		return s
	}
	ov := []comp{}
	for _, k := range overlap {
		ov = append(ov, mk(k))
	}
	rnd := mk(ASMLen)
	lo, hi := 32-12*sd, 32+12*sd
	const n = 8000
	h := (hi - lo) / n
	sum := 0.0
	for i := 0; i <= n; i++ {
		t := lo + float64(i)*h
		pdf := math.Exp(-0.5*((t-32)/sd)*((t-32)/sd)) / (sd * math.Sqrt(2*math.Pi))
		lg := float64(nRand) * math.Log(cdf(rnd, t))
		for _, c := range ov {
			lg += 2 * math.Log(cdf(c, t))
		}
		w := 2.0
		if i == 0 || i == n {
			w = 1
		} else if i%2 == 1 {
			w = 4
		}
		sum += w * pdf * math.Exp(lg)
	}
	return sum * h / 3
}
