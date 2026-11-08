// Astron Built this AP
// Project 022 | Rocket propulsion | Minimum-length planar nozzle by the method of characteristics
// Author: Aditya (aditya-astron)
// Date: 2026-11-08 IST
//
// Reference values: NACA Report 1135 tables (gamma = 1.4) and Anderson,
// Modern Compressible Flow, Example 11.1 (Me = 2.4 minimum-length nozzle).

package nozzle

import (
	"math"
	"testing"
)

const deg = math.Pi / 180

func near(t *testing.T, name string, got, want, tol float64) {
	t.Helper()
	if math.Abs(got-want) > tol {
		t.Errorf("%s = %.6f, want %.6f +/- %g", name, got, want, tol)
	}
}

func TestIsentropicTables(t *testing.T) {
	// NACA 1135, gamma = 1.4.
	near(t, "nu(2.0)", PrandtlMeyer(2.0, 1.4)/deg, 26.3798, 1e-4)
	near(t, "nu(2.4)", PrandtlMeyer(2.4, 1.4)/deg, 36.7465, 1e-4)
	near(t, "nu(3.0)", PrandtlMeyer(3.0, 1.4)/deg, 49.7573, 1e-4)
	near(t, "A/A*(2.0)", AreaRatio(2.0, 1.4), 1.6875, 1e-4)
	near(t, "A/A*(2.4)", AreaRatio(2.4, 1.4), 2.4031, 1e-4)
	near(t, "p/p0(2.0)", PressureRatio(2.0, 1.4), 0.12780, 1e-5)
	near(t, "A/A*(1)", AreaRatio(1, 1.22), 1, 1e-15)
	// Anderson Ex. 11.1: theta_max = nu(2.4)/2 = 18.375 deg (his table rounding).
	near(t, "theta_max", PrandtlMeyer(2.4, 1.4)/2/deg, 18.375, 2e-3)
}

func TestPrandtlMeyerInverse(t *testing.T) {
	for _, g := range []float64{1.13, 1.22, 1.4, 1.67} {
		for _, M := range []float64{1.0001, 1.05, 1.5, 2.4, 4, 8, 20} {
			back := MachFromPrandtlMeyer(PrandtlMeyer(M, g), g)
			if math.Abs(back-M) > 1e-9*M {
				t.Errorf("gamma %.2f: M %.4f -> %.12f", g, M, back)
			}
		}
	}
	if MachFromPrandtlMeyer(0, 1.4) != 1 {
		t.Error("nu = 0 must give M = 1")
	}
}

// Mass conservation: exit half-height converges to A/A*, monotonically and
// with error shrinking by more than 2x per grid doubling.
func TestExitAreaConverges(t *testing.T) {
	ar := AreaRatio(2.4, 1.4)
	prev := math.Inf(1)
	for _, n := range []int{15, 30, 60, 120, 240} {
		d, err := MinimumLengthNozzle(2.4, 1.4, n)
		if err != nil {
			t.Fatal(err)
		}
		e := math.Abs(d.ExitHeight()/ar - 1)
		if e > prev/2 {
			t.Errorf("n = %d: area error %.3e did not halve from %.3e", n, e, prev)
		}
		prev = e
	}
	if prev > 2e-4 {
		t.Errorf("area error at n = 240 is %.3e", prev)
	}
}

func TestExitStateAndWallShape(t *testing.T) {
	for _, c := range []struct{ Me, g float64 }{{2.4, 1.4}, {3.2, 1.22}, {4.0, 1.3}} {
		d, err := MinimumLengthNozzle(c.Me, c.g, 100)
		if err != nil {
			t.Fatal(err)
		}
		exit := d.Wall[len(d.Wall)-1]
		near(t, "exit theta", exit.Theta, 0, 1e-12)
		near(t, "exit M", exit.M, c.Me, 1e-9)
		near(t, "centreline exit M", d.Net[d.N-1][d.N-1].M, c.Me, 1e-9)
		// Corner turns to theta_max, then the wall straightens without
		// inflection: x and y increase, theta decreases.
		near(t, "corner angle", d.Wall[0].Theta, d.ThetaMax, 1e-12)
		for k := 1; k < len(d.Wall); k++ {
			a, b := d.Wall[k-1], d.Wall[k]
			if !(b.X > a.X && b.Y > a.Y && b.Theta < a.Theta) {
				t.Fatalf("Me %.1f: wall not monotone at %d", c.Me, k)
			}
		}
		// Every net point stays inside the nozzle.
		for i := range d.Net {
			for j := i; j < d.N; j++ {
				p := d.Net[i][j]
				if p.Y < -1e-12 || p.Y > d.ExitHeight() || p.X < 0 || p.X > d.Length() {
					t.Fatalf("net point (%d,%d) outside nozzle: %+v", i, j, p)
				}
			}
		}
	}
}

// Axial momentum: throat inflow plus wall pressure force equals exit flux.
// This only closes if the wall contour and wall pressures are both right.
func TestMomentumBalance(t *testing.T) {
	prev := math.Inf(1)
	for _, n := range []int{20, 40, 80, 160} {
		d, _ := MinimumLengthNozzle(3.2, 1.22, n)
		r := math.Abs(d.MomentumBalance())
		if r > prev/2 {
			t.Errorf("n = %d: residual %.3e did not halve from %.3e", n, r, prev)
		}
		prev = r
	}
	if prev > 5e-4 {
		t.Errorf("momentum residual at n = 160 is %.3e", prev)
	}
}

// For uniform axial exit flow, the momentum form of Cf equals the closed-form
// ideal Cf identically; and the MOC wall integral should approach it.
func TestThrustCoefficient(t *testing.T) {
	for _, c := range []struct{ Me, g, pa float64 }{{2.4, 1.4, 0}, {3.2, 1.22, 0}, {3.2, 1.22, 0.02}, {5, 1.15, 0.001}} {
		eps := AreaRatio(c.Me, c.g)
		pe := PressureRatio(c.Me, c.g)
		near(t, "Cf momentum vs closed form", ThrustCoefficient(c.Me, c.g, eps, c.pa), IdealCf(c.g, pe, eps, c.pa), 1e-12)
	}
	// Hand calc, gamma = 1.22, pe/p0 = 0.015236, eps = 8.304: Cf_vac = 1.708.
	d, _ := MinimumLengthNozzle(3.2, 1.22, 400)
	near(t, "MOC wall Cf_vac", d.WallThrustCoefficient(), 1.708, 1e-3)
	// Choked-throat sanity: Cf at eps = 1 (sonic exit, vacuum) is (1+g) p*/p0.
	g := 1.4
	near(t, "sonic Cf", ThrustCoefficient(1, g, 1, 0), (1+g)*PressureRatio(1, g), 1e-12)
}

func TestWedge(t *testing.T) {
	near(t, "lambda(15 deg)", WedgeDivergenceFactor(15*deg), 0.98862, 1e-5)
	near(t, "lambda -> 1", WedgeDivergenceFactor(1e-6), 1, 1e-12)
	// A wedge must underperform the shock-free uniform-exit nozzle.
	eps := AreaRatio(3.2, 1.22)
	if WedgeThrustCoefficient(3.2, 1.22, eps, 0, 15*deg) >= ThrustCoefficient(3.2, 1.22, eps, 0) {
		t.Error("wedge Cf should be below ideal")
	}
}

func TestBadInput(t *testing.T) {
	for _, c := range []struct {
		Me, g float64
		n     int
	}{{1, 1.4, 10}, {2, 1, 10}, {2, 1.4, 1}} {
		if _, err := MinimumLengthNozzle(c.Me, c.g, c.n); err == nil {
			t.Errorf("expected error for %+v", c)
		}
	}
}
