// Astron Built this AP
// Project 022 | Rocket propulsion | Minimum-length planar nozzle by the method of characteristics
// Author: Aditya (aditya-astron)
// Date: 2026-11-08 IST
//
// 1. Grid convergence on the textbook case (Me = 2.4, gamma = 1.4) against
//    the isentropic area ratio and an axial momentum balance.
// 2. A combustion-gas design (gamma = 1.22, Me = 3.2): contour, length and
//    thrust coefficient against a 15 deg planar wedge of the same area ratio.
//
// go run ./cmd/nozzle [-csv wall.csv]

package main

import (
	"flag"
	"fmt"
	"math"
	"os"

	nozzle "mocnozzle"
)

const deg = 180 / math.Pi

func must(d *nozzle.Design, err error) *nozzle.Design {
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
	return d
}

func main() {
	csvPath := flag.String("csv", "", "write the gamma=1.22 wall contour (x, y, theta_deg, M) to this file")
	flag.Parse()

	// Case 1: convergence.
	const me1, g1 = 2.4, 1.4
	ar := nozzle.AreaRatio(me1, g1)
	fmt.Printf("Case 1: Me = %.1f, gamma = %.2f\n", me1, g1)
	fmt.Printf("  nu(Me) = %.4f deg, theta_max = %.4f deg, isentropic A/A* = %.5f\n\n",
		nozzle.PrandtlMeyer(me1, g1)*deg, nozzle.PrandtlMeyer(me1, g1)/2*deg, ar)
	fmt.Println("     n    h_e/h*   area err    L/h*   momentum resid")
	var lengths []float64
	ns := []int{7, 15, 30, 60, 120, 240, 480}
	for _, n := range ns {
		d := must(nozzle.MinimumLengthNozzle(me1, g1, n))
		lengths = append(lengths, d.Length())
		fmt.Printf("  %4d  %8.5f  %+8.3f%%  %7.4f   %+.2e\n",
			n, d.ExitHeight(), 100*(d.ExitHeight()/ar-1), d.Length(), d.MomentumBalance())
	}
	// Observed order from the last three lengths, then Richardson.
	k := len(lengths)
	p := math.Log2((lengths[k-2] - lengths[k-3]) / (lengths[k-1] - lengths[k-2]))
	lInf := lengths[k-1] + (lengths[k-1]-lengths[k-2])/(math.Pow(2, p)-1)
	fmt.Printf("\n  observed order %.2f, extrapolated L/h* = %.4f\n\n", p, lInf)

	// Case 2: combustion-gas nozzle.
	const me2, g2, n2 = 3.2, 1.22, 200
	d := must(nozzle.MinimumLengthNozzle(me2, g2, n2))
	eps := nozzle.AreaRatio(me2, g2)
	pe := nozzle.PressureRatio(me2, g2)
	fmt.Printf("Case 2: Me = %.1f, gamma = %.2f, n = %d\n", me2, g2, n2)
	fmt.Printf("  theta_max      = %.3f deg\n", d.ThetaMax*deg)
	fmt.Printf("  eps (A/A*)     = %.4f isentropic, %.4f MOC\n", eps, d.ExitHeight())
	fmt.Printf("  pe/p0          = %.5f\n", pe)
	fmt.Printf("  length L/h*    = %.3f\n", d.Length())

	// Vacuum Cf from the MOC wall: throat momentum plus wall pressure integral.
	cfWall := d.WallThrustCoefficient()
	cfIdeal := nozzle.IdealCf(g2, pe, eps, 0)
	fmt.Printf("  Cf vac         = %.4f (throat + wall pressure), %.4f ideal 1-D\n", cfWall, cfIdeal)

	alpha := 15 / deg
	lam := nozzle.WedgeDivergenceFactor(alpha)
	cfWedge := nozzle.WedgeThrustCoefficient(me2, g2, eps, 0, alpha)
	fmt.Printf("\n  15 deg planar wedge, same eps:\n")
	fmt.Printf("    length L/h*  = %.3f (%.0f%% of MOC)\n", nozzle.WedgeLength(eps, alpha), 100*nozzle.WedgeLength(eps, alpha)/d.Length())
	fmt.Printf("    lambda       = %.4f\n", lam)
	fmt.Printf("    Cf vac       = %.4f (%.2f%% below MOC)\n", cfWedge, 100*(1-cfWedge/cfIdeal))

	if *csvPath != "" {
		f, err := os.Create(*csvPath)
		if err != nil {
			fmt.Fprintln(os.Stderr, err)
			os.Exit(1)
		}
		defer f.Close()
		fmt.Fprintln(f, "x_over_hstar,y_over_hstar,theta_deg,mach")
		for _, w := range d.Wall {
			fmt.Fprintf(f, "%.6f,%.6f,%.4f,%.5f\n", w.X, w.Y, w.Theta*deg, w.M)
		}
		fmt.Printf("\n  wall contour written to %s (%d points)\n", *csvPath, len(d.Wall))
	}
}
