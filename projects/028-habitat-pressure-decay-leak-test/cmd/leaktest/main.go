// Astron Built this AP
// Project 028 | Analog habitat systems | Habitat pressure-decay leak test
// Author: Aditya (aditya-astron)
// Date: 2026-11-10 IST
//
// Runs the reference 72 h leak test and compares three ways of reducing it:
// the raw pressure trace, the temperature-compensated two-point method over
// whole days, and the inventory fit of the full power law.
package main

import (
	"fmt"
	"math"
	"math/rand"

	"habitatleak"
)

func main() {
	sc := habitatleak.Reference()
	clean, mTrue := sc.TrueLog()
	log := sc.NoisyLog(clean, rand.New(rand.NewSource(28)))
	h := sc.Hab
	last := len(log.T) - 1

	const dpRef = 1000.0 // operating differential for the acceptance number
	rho := 101325.0 / (habitatleak.RAir * 293.15)
	qTrue := sc.Truth.Flow(dpRef)
	perDay := func(q float64) float64 { return q * rho * 86400 } // kg/day

	fmt.Println("Habitat pressure-decay leak test (reference scenario)")
	fmt.Printf("  volume %.0f m^3, sealed at +%.0f Pa, %.0f h log every %.0f s\n",
		h.Volume, sc.Dp0, sc.Duration/3600, sc.Step)
	fmt.Printf("  dp: %.0f Pa -> %.0f Pa, inventory %.2f kg -> %.2f kg\n",
		clean.Dp[0], clean.Dp[last], mTrue[0], mTrue[last])
	fmt.Printf("  temperature sensitivity p/T = %.0f Pa/K; barometer front -400 Pa\n\n", 101325/293.15)

	fmt.Printf("Truth: C = %.3e m^3/(s Pa^n), n = %.3f\n", sc.Truth.C, sc.Truth.N)
	fmt.Printf("  Q(%.0f Pa) = %.3e m^3/s = %.3f kg/day = %.3f %%/day of inventory\n\n",
		dpRef, qTrue, perDay(qTrue), 100*perDay(qTrue)/mTrue[0])

	// 1) Raw pressure trace, isothermal decay, n fixed at the true value so
	// the only error is the missing temperature and barometric correction.
	naive := habitatleak.NaiveDecay(log, h, sc.Truth.N)
	qN := naive.Flow(dpRef)
	fmt.Println("1) Pressure trace only (no T or barometer correction, n given):")
	fmt.Printf("   Q(%.0f Pa) = %.3e m^3/s  error %+.1f %%\n\n", dpRef, qN, 100*(qN/qTrue-1))

	// 2) Two-point inventory method over each whole day, so the diurnal
	// temperature swing mostly cancels.
	fmt.Println("2) Two-point temperature-compensated method, 24 h windows:")
	per := int(86400 / sc.Step)
	for d := 0; d < 3; d++ {
		a, b := d*per, (d+1)*per
		mdot, dpm := habitatleak.EndpointLeak(log, h, a, b)
		mdotTrue := (mTrue[a] - mTrue[b]) / (log.T[b] - log.T[a])
		fmt.Printf("   day %d: mdot = %.4f kg/day at mean dp %4.0f Pa  (true mass loss %.4f, error %+.1f %%)\n",
			d+1, mdot*86400, dpm, mdotTrue*86400, 100*(mdot/mdotTrue-1))
	}
	fmt.Println("   (one rate per window: no exponent, so no extrapolation to other dp)")
	fmt.Println()

	// 3) Inventory fit of C and n with the measured weather as input.
	fit, err := h.FitLeak(log, habitatleak.LeakLaw{C: 1e-6, N: 0.75}, 1)
	if err != nil {
		fmt.Println("fit failed:", err)
		return
	}
	q, sq := fit.FlowAt(dpRef)
	sn := math.Sqrt(fit.Cov[1][1])
	fmt.Printf("3) Inventory fit (Levenberg-Marquardt, %d iterations, residual %.1f g):\n", fit.Iter, 1000*fit.Sigma)
	fmt.Printf("   n = %.3f +/- %.3f   (truth %.3f)\n", fit.Law.N, sn, sc.Truth.N)
	fmt.Printf("   C = %.3e  (truth %.3e)\n", fit.Law.C, sc.Truth.C)
	fmt.Printf("   Q(%.0f Pa) = %.4e +/- %.1e m^3/s  error %+.2f %% (%.1f sigma)\n",
		dpRef, q, sq, 100*(q/qTrue-1), (q-qTrue)/sq)
	fmt.Printf("   leak = %.3f kg/day = %.3f %%/day at %.0f Pa\n", perDay(q), 100*perDay(q)/fit.M0, dpRef)

	for _, dp := range []float64{dpRef, sc.Dp0} {
		qa := fit.Law.Flow(dp)
		a := habitatleak.EquivalentArea(qa, dp, rho, 0.61)
		fmt.Printf("   equivalent orifice (Cd 0.61) at %4.0f Pa: %.3f mm^2, d = %.3f mm\n",
			dp, a*1e6, 2*math.Sqrt(a/math.Pi)*1e3)
	}
}
