// Astron Built this AP
// Project 045 | Thermal analysis | Multilayer insulation layer-by-layer heat leak
// Author: Aditya (aditya-astron)
// Date: 2026-11-17 IST

package main

import (
	"fmt"
	"math"
)

func must(s *Solution, err error) *Solution {
	if err != nil {
		panic(err)
	}
	return s
}

func maxDiff(a, b []float64) float64 {
	d := 0.0
	for i := range a {
		d = math.Max(d, math.Abs(a[i]-b[i]))
	}
	return d
}

// flightBlanket: 18 DAM shields plus an outer cover, inner layer held at the
// 293 K bus wall. Cover optics are assumed values for an aged Kapton-faced
// outer layer (inner face aluminised).
func flightBlanket(flux float64) *Blanket {
	b := UniformBlanket(18, 0.03, 15, 1e-7, 293, 0)
	b.Radiating = true
	b.Alpha, b.EpsOut, b.Flux = 0.40, 0.75, flux
	return b
}

func main() {
	fmt.Println("== 1. Pure radiation, 20 shields, eps 0.03, 293 K -> 150 K ==")
	pr := UniformBlanket(20, 0.03, 0, 0, 293, 150)
	sr := must(pr.SolveNewton())
	qHand := sigma * (pow4(293) - pow4(150)) / (21 * (2/0.03 - 1))
	fmt.Printf("Newton q = %.6f W/m^2 (%d iters), hand sigma dT4/((n+1)(2/eps-1)) = %.6f, rel err %.1e\n",
		sr.Q, sr.Iters, qHand, math.Abs(sr.Q-qHand)/qHand)
	fmt.Printf("mid-stack shield T10 = %.3f K, T^4-linear prediction %.3f K\n",
		sr.T[10], math.Pow(pow4(293)-10.0/21*(pow4(293)-pow4(150)), 0.25))

	fmt.Println("\n== 2. Lab boundary 300 K -> 77 K, 20 shields, Nbar 20 /cm, 1e-6 torr ==")
	lab := UniformBlanket(20, 0.031, 20, 1e-6, 300, 77)
	sn := must(lab.SolveNewton())
	cf := lab.ClosedForm()
	sh := lab.SolveShooting()
	fmt.Printf("Newton      q = %.6f W/m^2 (%d iters)\n", sn.Q, sn.Iters)
	fmt.Printf("closed form q = %.6f  (F(Th)-F(Tc))/(n+1), max |dT| %.2e K\n", cf.Q, maxDiff(sn.T, cf.T))
	fmt.Printf("shooting    q = %.6f  max |dT| %.2e K\n", sh.Q, maxDiff(sn.T, sh.T))
	r := lab.radFactor(0) * (pow4(300) - pow4(77)) / 21
	s := 0.5 * lab.hs() * (300*300 - 77*77) / 21
	g := sn.Q - r - s
	lr, ls, lg := Lockheed(20, 0.031, 20, 1e-6, 300, 77)
	fmt.Printf("model    : rad %.4f  solid %.4f  gas %.4f  total %.4f W/m^2  eps* %.5f\n",
		r, s, g, sn.Q, EffectiveEmittance(sn.Q, 300, 77))
	fmt.Printf("Lockheed : rad %.4f  solid %.4f  gas %.4f  total %.4f W/m^2  eps* %.5f\n",
		lr, ls, lg, lr+ls+lg, EffectiveEmittance(lr+ls+lg, 300, 77))
	fmt.Println("node   T [K]")
	for _, i := range []int{0, 1, 5, 10, 15, 20, 21} {
		fmt.Printf("%4d  %7.2f\n", i, sn.T[i])
	}

	fmt.Println("\n== 3. Non-uniform stack: blanket on 0.85 paint, 3 crinkled inner shields eps 0.08 ==")
	nu := UniformBlanket(15, 0.03, 18, 1e-6, 293, 120)
	nu.EpsCold[0] = 0.85
	for i := 1; i <= 3; i++ {
		nu.EpsHot[i], nu.EpsCold[i] = 0.08, 0.08
	}
	nn := must(nu.SolveNewton())
	ns := nu.SolveShooting()
	fmt.Printf("Newton q = %.6f (%d iters), shooting q = %.6f, rel %.1e, max |dT| %.2e K\n",
		nn.Q, nn.Iters, ns.Q, math.Abs(nn.Q-ns.Q)/nn.Q, maxDiff(nn.T, ns.T))
	fmt.Printf("gap fluxes: first %.6f last %.6f\n", nn.Q, func() float64 {
		q, _, _ := nu.GapFlux(15, nn.T[15], nn.T[16])
		return q
	}())

	fmt.Println("\n== 4. Flight blanket on a 293 K wall: 18 shields, cover alpha 0.40 eps 0.75 ==")
	for _, c := range []struct {
		name string
		s    float64
	}{{"eclipse", 0}, {"sun 45 deg", 1361 * math.Cos(math.Pi/4)}, {"sun normal", 1361}} {
		b := flightBlanket(c.s)
		x := must(b.SolveNewton())
		y := b.SolveShooting()
		fmt.Printf("%-11s cover %6.1f K  q into blanket %+8.4f W/m^2  (shoot %+8.4f)  eps* vs cover %.5f\n",
			c.name, x.T[19], x.Q, y.Q, EffectiveEmittance(x.Q, 293, x.T[19]))
	}
	// Normal solar flux at which the blanket carries no heat either way.
	lo, hi := 0.0, 1361.0
	for i := 0; i < 100; i++ {
		mid := 0.5 * (lo + hi)
		if must(flightBlanket(mid).SolveNewton()).Q > 0 {
			lo = mid
		} else {
			hi = mid
		}
	}
	zero := 0.5 * (lo + hi)
	fmt.Printf("zero-leak absorbed-sun flux %.2f W/m^2 (cover at %.2f K, radiative check %.2f K)\n",
		zero, must(flightBlanket(zero).SolveNewton()).T[19], math.Pow(0.40*zero/(0.75*sigma), 0.25))
	fmt.Printf("eclipse heater make-up for 2.5 m^2 of blanket: %.2f W\n",
		2.5*must(flightBlanket(0).SolveNewton()).Q)

	fmt.Println("\n== 5. Fixed 1.0 cm blanket, 293 K -> 150 K: how many shields? ==")
	f := func(n float64) float64 { return FixedThicknessFlux(n, 1.0, 0.03, 1e-6, 293, 150) }
	nStar := GoldenMin(f, 1, 200)
	best := 1
	for n := 1; n <= 200; n++ {
		if f(float64(n)) < f(float64(best)) {
			best = n
		}
	}
	fmt.Printf("continuous optimum n* = %.2f (q %.4f), integer brute force n = %d (q %.4f W/m^2)\n",
		nStar, f(nStar), best, f(float64(best)))
	fmt.Println("   n   Nbar/cm   q W/m^2   layer-model Newton")
	for _, n := range []int{5, 10, 20, best, 60, 100} {
		b := UniformBlanket(n, 0.03, float64(n), 1e-6, 293, 150)
		fmt.Printf("%4d   %6.1f   %8.4f   %8.4f\n", n, float64(n), f(float64(n)), must(b.SolveNewton()).Q)
	}
}
