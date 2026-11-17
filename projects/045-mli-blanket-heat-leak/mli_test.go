// Astron Built this AP
// Project 045 | Thermal analysis | Multilayer insulation layer-by-layer heat leak
// Author: Aditya (aditya-astron)
// Date: 2026-11-17 IST

package main

import (
	"math"
	"testing"
)

func relErr(a, b float64) float64 { return math.Abs(a-b) / math.Abs(b) }

// N grey shields between grey plates: q = sigma dT^4 / ((n+1)(2/eps-1)), and
// T^4 falls linearly through the stack.
func TestPureRadiationHandFormula(t *testing.T) {
	for _, n := range []int{1, 5, 20, 60} {
		b := UniformBlanket(n, 0.05, 0, 0, 293, 150)
		s, err := b.SolveNewton()
		if err != nil {
			t.Fatal(err)
		}
		want := sigma * (pow4(293) - pow4(150)) / (float64(n+1) * (2/0.05 - 1))
		if e := relErr(s.Q, want); e > 1e-12 {
			t.Errorf("n=%d q=%g want %g rel %g", n, s.Q, want, e)
		}
		for k, tk := range s.T {
			t4 := pow4(293) - float64(k)/float64(n+1)*(pow4(293)-pow4(150))
			if math.Abs(pow4(tk)-t4)/t4 > 1e-11 {
				t.Errorf("n=%d node %d T^4 not linear", n, k)
			}
		}
	}
}

// With radiation off, the per-gap spacer conductance must reproduce the
// Lockheed solid term exactly (the telescoping argument in mli.go).
func TestSolidConductionMatchesLockheed(t *testing.T) {
	b := UniformBlanket(25, 0.03, 22, 0, 300, 77)
	b.NoRad = true
	s, err := b.SolveNewton()
	if err != nil {
		t.Fatal(err)
	}
	_, qs, _ := Lockheed(25, 0.03, 22, 0, 300, 77)
	if e := relErr(s.Q, qs); e > 1e-12 {
		t.Errorf("solid q %g Lockheed %g rel %g", s.Q, qs, e)
	}
}

func TestNewtonMatchesClosedFormAndShooting(t *testing.T) {
	b := UniformBlanket(20, 0.031, 20, 1e-6, 300, 77)
	sn, err := b.SolveNewton()
	if err != nil {
		t.Fatal(err)
	}
	if sn.Iters > 10 {
		t.Errorf("Newton took %d iterations", sn.Iters)
	}
	cf := b.ClosedForm()
	sh := b.SolveShooting()
	if relErr(sn.Q, cf.Q) > 1e-12 || maxDiff(sn.T, cf.T) > 1e-9 {
		t.Errorf("closed form mismatch q %g vs %g", sn.Q, cf.Q)
	}
	if relErr(sn.Q, sh.Q) > 1e-10 || maxDiff(sn.T, sh.T) > 1e-8 {
		t.Errorf("shooting mismatch q %g vs %g", sn.Q, sh.Q)
	}
	// Lab-grade DAM blankets measure eps* of order 1e-3; the ideal-shield
	// model should sit within 15 % of the Lockheed fit at its own conditions.
	lr, ls, lg := Lockheed(20, 0.031, 20, 1e-6, 300, 77)
	if e := relErr(sn.Q, lr+ls+lg); e > 0.15 {
		t.Errorf("model vs Lockheed rel %g", e)
	}
}

func TestNonUniformStackFluxBalance(t *testing.T) {
	b := UniformBlanket(15, 0.03, 18, 1e-6, 293, 120)
	b.EpsCold[0] = 0.85
	for i := 1; i <= 3; i++ {
		b.EpsHot[i], b.EpsCold[i] = 0.08, 0.08
	}
	sn, err := b.SolveNewton()
	if err != nil {
		t.Fatal(err)
	}
	for k := 0; k <= b.N; k++ {
		q, _, _ := b.GapFlux(k, sn.T[k], sn.T[k+1])
		if relErr(q, sn.Q) > 1e-11 {
			t.Errorf("gap %d flux %g != %g", k, q, sn.Q)
		}
	}
	sh := b.SolveShooting()
	if relErr(sn.Q, sh.Q) > 1e-10 {
		t.Errorf("shooting q %g vs Newton %g", sh.Q, sn.Q)
	}
	// Degrading layers can only raise the leak relative to the clean stack.
	clean, _ := UniformBlanket(15, 0.03, 18, 1e-6, 293, 120).SolveNewton()
	if sn.Q <= clean.Q {
		t.Errorf("degraded stack %g not above clean %g", sn.Q, clean.Q)
	}
}

func TestRadiatingCover(t *testing.T) {
	for _, s := range []float64{0, 500, 1361} {
		b := flightBlanket(s)
		sn, err := b.SolveNewton()
		if err != nil {
			t.Fatal(err)
		}
		sh := b.SolveShooting()
		if math.Abs(sn.Q-sh.Q) > 1e-10 {
			t.Errorf("S=%g Newton %g shooting %g", s, sn.Q, sh.Q)
		}
		// Cover energy balance: blanket flux + absorbed sun = emitted.
		te := sn.T[b.N+1]
		if math.Abs(sn.Q+b.Alpha*s-b.EpsOut*sigma*pow4(te)) > 1e-9 {
			t.Errorf("S=%g cover balance off", s)
		}
	}
	// Zero leak exactly when the cover sits at wall temperature.
	s0 := 0.75 * sigma * pow4(293) / 0.40
	sn, _ := flightBlanket(s0).SolveNewton()
	if math.Abs(sn.Q) > 1e-9 || math.Abs(sn.T[19]-293) > 1e-7 {
		t.Errorf("zero-leak case q=%g Tcover=%g", sn.Q, sn.T[19])
	}
}

func TestFixedThicknessOptimum(t *testing.T) {
	f := func(n float64) float64 { return FixedThicknessFlux(n, 1.0, 0.03, 1e-6, 293, 150) }
	nStar := GoldenMin(f, 1, 200)
	best := 1
	for n := 1; n <= 200; n++ {
		if f(float64(n)) < f(float64(best)) {
			best = n
		}
	}
	if best != int(math.Floor(nStar)) && best != int(math.Ceil(nStar)) {
		t.Errorf("integer optimum %d not adjacent to %g", best, nStar)
	}
	// Stationarity of the continuous optimum by central difference.
	h := 1e-3
	if d := (f(nStar+h) - f(nStar-h)) / (2 * h); math.Abs(d) > 1e-6 {
		t.Errorf("df/dn at optimum = %g", d)
	}
	// Closed form must agree with the Newton layer model at the integer optimum.
	s, _ := UniformBlanket(best, 0.03, float64(best), 1e-6, 293, 150).SolveNewton()
	if relErr(s.Q, f(float64(best))) > 1e-12 {
		t.Errorf("layer model %g vs closed form %g", s.Q, f(float64(best)))
	}
}
