// Astron Built this AP
// Project 028 | Analog habitat systems | Habitat pressure-decay leak test
// Author: Aditya (aditya-astron)
// Date: 2026-11-10 IST
//
// Verification: hand calculations, the ODE integrator against closed-form and
// quadrature decay times, exact inventory compensation of a sealed hull, RK4
// order, noiseless parameter recovery, and a Monte Carlo check that the fit
// is unbiased and its reported covariance matches the actual scatter.
package habitatleak

import (
	"math"
	"math/rand"
	"testing"
)

func relErr(a, b float64) float64 { return math.Abs(a/b - 1) }

type still struct{ pe, ta float64 }

func (s still) At(float64) (float64, float64) { return s.pe, s.ta }

// Hand calculation: Q = 1e-5 m^3/s at 100 Pa, rho 1.2 kg/m^3, Cd 0.61:
// v = sqrt(2*100/1.2) = 12.90994 m/s, A = 1e-5/(0.61*12.90994) = 1.26983e-6 m^2.
func TestEquivalentAreaHand(t *testing.T) {
	a := EquivalentArea(1e-5, 100, 1.2, 0.61)
	if relErr(a, 1.26983e-6) > 1e-5 {
		t.Fatalf("A = %.6e, want 1.26983e-6", a)
	}
	law := LeakLaw{C: 2e-7, N: 0.65}
	// 2e-7 * 1000^0.65 = 2e-7 * 89.1251 = 1.78250e-5
	if relErr(law.Flow(1000), 1.78250e-5) > 1e-5 || law.Flow(-1000) != -law.Flow(1000) {
		t.Fatalf("Flow(1000) = %.6e", law.Flow(1000))
	}
}

// A sealed hull (C = 0) sees hundreds of Pa of weather-driven dp swing, yet
// the inventory computed from the log must stay constant to round-off.
// This is the whole premise of the method.
func TestSealedHullInventoryConstant(t *testing.T) {
	sc := Reference()
	sc.Truth = LeakLaw{C: 0, N: 0.62}
	s, _ := sc.TrueLog()
	m := sc.Hab.Mass(s)
	lo, hi := s.Dp[0], s.Dp[0]
	for k := range m {
		if relErr(m[k], m[0]) > 1e-13 {
			t.Fatalf("inventory drifted at k=%d: %.15g vs %.15g", k, m[k], m[0])
		}
		lo, hi = math.Min(lo, s.Dp[k]), math.Max(hi, s.Dp[k])
	}
	// 2.6 K peak-to-peak air temperature is ~900 Pa, the front adds 400 Pa.
	if hi-lo < 800 {
		t.Fatalf("dp swing only %.0f Pa; the scenario should be weather-dominated", hi-lo)
	}
	t.Logf("sealed hull: dp swings %.0f Pa with zero leak, inventory constant", hi-lo)
}

// With constant weather, dp' = -(p_e + dp) C dp^n / V, so the time to fall
// from dp0 to dp1 is t = int_{dp1}^{dp0} V / (C (p_e + x) x^n) dx.
func decayTime(V, C, n, pe, dp0, dp1 float64) float64 {
	const N = 4000 // composite Simpson, even
	h := (dp0 - dp1) / N
	f := func(x float64) float64 { return V / (C * (pe + x) * math.Pow(x, n)) }
	s := f(dp1) + f(dp0)
	for i := 1; i < N; i++ {
		w := 2.0
		if i%2 == 1 {
			w = 4
		}
		s += w * f(dp1+float64(i)*h)
	}
	return s * h / 3
}

func dpAfter(h Habitat, law LeakLaw, w still, dp0, tEnd float64, sub int) float64 {
	m0 := (w.pe + dp0) * h.Volume / (RAir * w.ta)
	m := h.Integrate(law, w, m0, []float64{0, tEnd}, sub)
	return m[1]*RAir*w.ta/h.Volume - w.pe
}

func TestIntegratorAgainstQuadrature(t *testing.T) {
	h := Habitat{Volume: 300, TExt: 288.15}
	w := still{pe: 101325, ta: 293.15}
	law := LeakLaw{C: 2.3e-7, N: 0.62}
	tEnd := 3 * day
	dp1 := dpAfter(h, law, w, 3000, tEnd, 3*864)
	tq := decayTime(h.Volume, law.C, law.N, w.pe, 3000, dp1)
	if relErr(tq, tEnd) > 1e-8 {
		t.Fatalf("quadrature time %.6f s vs integrated %.6f s", tq, tEnd)
	}
	t.Logf("n=0.62: dp(72 h) = %.4f Pa, quadrature time error %.2e", dp1, relErr(tq, tEnd))
}

// n = 1 has a closed form:
// t = V/(C p_e) [ ln(dp0/dp1) - ln((p_e+dp0)/(p_e+dp1)) ].
func TestIntegratorLaminarClosedForm(t *testing.T) {
	h := Habitat{Volume: 300, TExt: 288.15}
	w := still{pe: 101325, ta: 293.15}
	law := LeakLaw{C: 5e-9, N: 1}
	tEnd := 2 * day
	dp1 := dpAfter(h, law, w, 3000, tEnd, 2*864)
	tc := h.Volume / (law.C * w.pe) * (math.Log(3000/dp1) - math.Log((w.pe+3000)/(w.pe+dp1)))
	if relErr(tc, tEnd) > 1e-9 {
		t.Fatalf("closed-form time %.6f s vs %.6f s", tc, tEnd)
	}
}

// Fourth order: halving the step on the time-varying weather cuts the error
// by about 16.
func TestRK4Order(t *testing.T) {
	sc := Reference()
	times := []float64{0, sc.Duration}
	pe0, ta0 := Weather{}.At(0)
	m0 := (pe0 + sc.Dp0) * sc.Hab.Volume / (RAir * ta0)
	ref := sc.Hab.Integrate(sc.Truth, Weather{}, m0, times, 20000)[1]
	e1 := math.Abs(sc.Hab.Integrate(sc.Truth, Weather{}, m0, times, 60)[1] - ref)
	e2 := math.Abs(sc.Hab.Integrate(sc.Truth, Weather{}, m0, times, 120)[1] - ref)
	ratio := e1 / e2
	if ratio < 12 || ratio > 20 {
		t.Fatalf("error ratio %.2f, want ~16 (e1 %.3e, e2 %.3e)", ratio, e1, e2)
	}
	t.Logf("RK4 error ratio on step halving: %.2f", ratio)
}

func TestInflowUsesExteriorDensity(t *testing.T) {
	h := Habitat{Volume: 300, TExt: 288.15}
	law := LeakLaw{C: 2.3e-7, N: 0.62}
	pe, ta := 101325.0, 293.15
	m := (pe - 500) * h.Volume / (RAir * ta) // 500 Pa underpressure
	got := h.MassRate(law, m, pe, ta)
	want := pe / (RAir * h.TExt) * law.C * math.Pow(500, 0.62)
	if relErr(got, want) > 1e-12 {
		t.Fatalf("inflow rate %.6e, want %.6e", got, want)
	}
}

// Noiseless log: the only model error is linear interpolation of the logged
// weather, so C and n must come back almost exactly.
func TestFitNoiseless(t *testing.T) {
	sc := Reference()
	s, _ := sc.TrueLog()
	fit, err := sc.Hab.FitLeak(s, LeakLaw{C: 1e-6, N: 0.9}, 4)
	if err != nil {
		t.Fatal(err)
	}
	if relErr(fit.Law.N, sc.Truth.N) > 1e-3 || relErr(fit.Law.C, sc.Truth.C) > 5e-3 {
		t.Fatalf("recovered C=%.6e n=%.6f, truth C=%.6e n=%.6f",
			fit.Law.C, fit.Law.N, sc.Truth.C, sc.Truth.N)
	}
	q, _ := fit.FlowAt(1000)
	if relErr(q, sc.Truth.Flow(1000)) > 1e-4 {
		t.Fatalf("Q(1000) off by %.2e", relErr(q, sc.Truth.Flow(1000)))
	}
	t.Logf("noiseless: n err %.1e, C err %.1e, Q(1000) err %.1e, %d iterations",
		relErr(fit.Law.N, sc.Truth.N), relErr(fit.Law.C, sc.Truth.C),
		relErr(q, sc.Truth.Flow(1000)), fit.Iter)
}

// The uncorrected pressure trace is biased by the weather; the inventory fit
// is not. Same log for both.
func TestNaiveIsBiasedFitIsNot(t *testing.T) {
	sc := Reference()
	clean, _ := sc.TrueLog()
	s := sc.NoisyLog(clean, rand.New(rand.NewSource(1)))
	qt := sc.Truth.Flow(1000)
	naive := NaiveDecay(s, sc.Hab, sc.Truth.N).Flow(1000)
	fit, err := sc.Hab.FitLeak(s, LeakLaw{C: 1e-6, N: 0.75}, 1)
	if err != nil {
		t.Fatal(err)
	}
	q, _ := fit.FlowAt(1000)
	if relErr(naive, qt) < 0.05 || relErr(q, qt) > 0.01 {
		t.Fatalf("naive err %.3f, fit err %.4f", relErr(naive, qt), relErr(q, qt))
	}
	t.Logf("Q(1000 Pa) error: pressure-only %.1f %%, inventory fit %.2f %%",
		100*relErr(naive, qt), 100*relErr(q, qt))
}

// Monte Carlo over instrument noise: the fit must be unbiased and the
// covariance it reports must describe the real scatter.
func TestMonteCarloCovariance(t *testing.T) {
	if testing.Short() {
		t.Skip("Monte Carlo skipped in -short mode")
	}
	sc := Reference()
	clean, _ := sc.TrueLog()
	rng := rand.New(rand.NewSource(2026))
	const runs = 150
	qt := sc.Truth.Flow(1000)
	var sumQ, sumQ2, sumSig, sumN, sumN2, sumSN float64
	inside := 0
	for r := 0; r < runs; r++ {
		s := sc.NoisyLog(clean, rng)
		fit, err := sc.Hab.FitLeak(s, LeakLaw{C: 1e-6, N: 0.75}, 1)
		if err != nil {
			t.Fatal(err)
		}
		q, sq := fit.FlowAt(1000)
		e := q/qt - 1
		sumQ += e
		sumQ2 += e * e
		sumSig += sq / qt
		dn := fit.Law.N - sc.Truth.N
		sn := math.Sqrt(fit.Cov[1][1])
		sumN += dn
		sumN2 += dn * dn
		sumSN += sn
		if math.Abs(dn) < 2*sn {
			inside++
		}
	}
	mean := sumQ / runs
	sd := math.Sqrt(sumQ2/runs - mean*mean)
	pred := sumSig / runs
	meanN := sumN / runs
	sdN := math.Sqrt(sumN2/runs - meanN*meanN)
	predN := sumSN / runs
	t.Logf("Q(1000): bias %+.4f %%, empirical sd %.4f %%, predicted %.4f %%", 100*mean, 100*sd, 100*pred)
	t.Logf("n:       bias %+.5f, empirical sd %.5f, predicted %.5f, %d/%d inside 2 sigma",
		meanN, sdN, predN, inside, runs)

	if math.Abs(mean) > 4*sd/math.Sqrt(runs) {
		t.Errorf("Q bias %.4f %% exceeds 4 standard errors", 100*mean)
	}
	if math.Abs(meanN) > 4*sdN/math.Sqrt(runs) {
		t.Errorf("n bias %.5f exceeds 4 standard errors", meanN)
	}
	if sd/pred < 0.75 || sd/pred > 1.33 || sdN/predN < 0.75 || sdN/predN > 1.33 {
		t.Errorf("covariance mismatch: Q sd/pred %.2f, n sd/pred %.2f", sd/pred, sdN/predN)
	}
	if inside < int(0.88*runs) {
		t.Errorf("only %d/%d runs inside 2 sigma for n", inside, runs)
	}
}
