// Astron Built this AP
// Project 059 | Telemetry framing | Optimal ASM frame synchronization: Massey rule vs correlation
// Author: Aditya (aditya-astron)
// Date: 2026-11-23 IST

// Study driver: ASM overlap sidelobes, hard-threshold detector design, and
// one-shot acquisition of a CCSDS CADU at the low channel SNRs where LDPC
// and turbo codes operate, comparing the sync rules by Monte Carlo against
// the exact-overlap independence approximation.

package main

import (
	"flag"
	"fmt"
	"math"
	"math/rand/v2"
	"sync"
)

// CADU of an RS(255,223) I=5 frame: 32-bit ASM + 1275 octets.
const caduBits = 32 + 1275*8

func lin(db float64) float64 { return math.Pow(10, db/10) }

func sidelobes() {
	fmt.Println("== ASM 0x1ACFFC1D overlap sidelobes ==")
	fmt.Println("shift k  overlap  mismatches m_k  overlap corr (32-k-2m_k)")
	worst := 0
	for k := 1; k < ASMLen; k++ {
		m := OverlapMismatch(k)
		c := ASMLen - k - 2*m
		if c > worst {
			worst = c
		}
		if k <= 8 || k%4 == 0 {
			fmt.Printf("%7d  %7d  %14d  %+10d\n", k, ASMLen-k, m, c)
		}
	}
	fmt.Printf("largest aperiodic sidelobe over k = 1..31: %+d (main peak 32)\n\n", worst)
}

func thresholdTable() {
	fmt.Println("== Hard-decision threshold detector: declare sync if d_H <= E ==")
	fmt.Println("P_F random = sum_{i<=E} C(32,i)/2^32; worst overlap = max_k P(D_k <= E), noiseless")
	nRand := caduBits - 1 - 2*(ASMLen-1)
	for _, esdb := range []float64{-2, 2, 6} {
		p := BER(lin(esdb))
		fmt.Printf("\nEs/N0 = %+.0f dB, BER = %.3e\n", esdb, p)
		fmt.Println(" E   P_F random    worst overlap  false/frame    P_miss")
		for E := 0; E <= 6; E++ {
			pf := 0.0
			for i := 0; i <= E; i++ {
				pf += Binom(ASMLen, i)
			}
			pf /= math.Pow(2, ASMLen)
			worstOv, sumOv := 0.0, 0.0
			for k := 1; k < ASMLen; k++ {
				pmf := DistancePMF(k, p)
				q := 0.0
				for i := 0; i <= E && i < len(pmf); i++ {
					q += pmf[i]
				}
				sumOv += 2 * q
				pn := DistancePMF(k, 0)
				qn := 0.0
				for i := 0; i <= E && i < len(pn); i++ {
					qn += pn[i]
				}
				worstOv = math.Max(worstOv, qn)
			}
			pt := DistancePMF(0, p)
			miss := 0.0
			for i := E + 1; i <= ASMLen; i++ {
				miss += pt[i]
			}
			fmt.Printf("%2d  %11.3e  %13.3e  %11.3e  %9.3e\n",
				E, pf, worstOv, float64(nRand)*pf+sumOv, miss)
		}
	}
	fmt.Println()
}

func monteCarlo(L int, esn0 float64, trials int, seed uint64) [NumRules]int {
	const workers = 4
	var mu sync.Mutex
	var fails [NumRules]int
	var wg sync.WaitGroup
	for w := 0; w < workers; w++ {
		wg.Add(1)
		go func(w int) {
			defer wg.Done()
			rng := rand.New(rand.NewPCG(seed, uint64(w)+1))
			s := NewSearcher(L, esn0)
			var local [NumRules]int
			for t := w; t < trials; t += workers {
				ok := s.Trial(rng)
				for r := range ok {
					if !ok[r] {
						local[r]++
					}
				}
			}
			mu.Lock()
			for r := range local {
				fails[r] += local[r]
			}
			mu.Unlock()
		}(w)
	}
	wg.Wait()
	return fails
}

// crossing returns the Es/N0 (dB) at which a decreasing failure curve,
// interpolated in log10(P), reaches target.
func crossing(db, pf []float64, target float64) float64 {
	for i := 0; i+1 < len(db); i++ {
		a, b := pf[i], pf[i+1]
		if a >= target && b <= target && a > 0 && b > 0 {
			f := (math.Log10(a) - math.Log10(target)) / (math.Log10(a) - math.Log10(b))
			return db[i] + f*(db[i+1]-db[i])
		}
	}
	return math.NaN()
}

func acquisition(trials int) {
	L := caduBits
	fmt.Printf("== One-shot acquisition over one CADU (L = %d candidate offsets) ==\n", L)
	fmt.Printf("P(fail) = P(argmax metric != true offset); %d trials per point\n", trials)
	fmt.Println("analytic: independent windows with exact overlap pmfs (hard, corr)")
	snrs := []float64{-8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4}
	curves := make([][]float64, NumRules)
	fmt.Printf("%7s  %9s %9s", "Es/N0", "hard th", "corr th")
	for r := Rule(0); r < NumRules; r++ {
		fmt.Printf(" %16s", RuleNames[r])
	}
	fmt.Println()
	for i, db := range snrs {
		esn0 := lin(db)
		fails := monteCarlo(L, esn0, trials, 0x5907+uint64(i))
		ph := 1 - HardAcqAnalytic(L, BER(esn0))
		pc := 1 - CorrAcqAnalytic(L, esn0)
		fmt.Printf("%+5.0f dB  %9.4f %9.4f", db, ph, pc)
		for r := range fails {
			pf := float64(fails[r]) / float64(trials)
			curves[r] = append(curves[r], pf)
			fmt.Printf(" %16.4f", pf)
		}
		fmt.Println()
	}
	fmt.Println()
	for _, target := range []float64{0.1, 0.01} {
		ref := crossing(snrs, curves[RuleCorr], target)
		fmt.Printf("Es/N0 for P(fail) = %.2f (interpolated):\n", target)
		for r := Rule(0); r < NumRules; r++ {
			x := crossing(snrs, curves[r], target)
			fmt.Printf("  %-17s %+6.2f dB   (%+.2f dB vs corr)\n", RuleNames[r], x, x-ref)
		}
	}
	fmt.Println()
}

func main() {
	trials := flag.Int("trials", 5000, "Monte Carlo trials per SNR point")
	flag.Parse()
	sidelobes()
	thresholdTable()
	acquisition(*trials)
}
