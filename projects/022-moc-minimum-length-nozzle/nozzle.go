// Astron Built this AP
// Project 022 | Rocket propulsion | Minimum-length planar nozzle by the method of characteristics
// Author: Aditya (aditya-astron)
// Date: 2026-11-08 IST

// Package nozzle designs the shortest shock-free planar supersonic nozzle that
// delivers uniform, axial flow at a chosen exit Mach number.
//
// The throat is a sharp corner. A centred Prandtl-Meyer fan turns the flow at
// the corner up to theta_max = nu(Me)/2; the fan waves reflect off the
// centreline, and the wall downstream is shaped to cancel each reflected wave
// as it arrives. In planar, homentropic, irrotational flow the Riemann
// invariants
//
//	K- = theta + nu   (constant along C- characteristics, slope theta - mu)
//	K+ = theta - nu   (constant along C+ characteristics, slope theta + mu)
//
// make the flow state at every grid point algebraic; only the geometry needs
// the unit processes. Throat half-height is 1, so wall y is also h/h*.
package nozzle

import (
	"errors"
	"math"
)

// PrandtlMeyer returns nu(M) in radians for M >= 1.
func PrandtlMeyer(M, g float64) float64 {
	a := math.Sqrt((g + 1) / (g - 1))
	b := math.Sqrt(M*M - 1)
	return a*math.Atan(b/a) - math.Atan(b)
}

// MachFromPrandtlMeyer inverts nu(M). nu is monotone in M, so Newton with a
// bisection fallback converges from any start.
func MachFromPrandtlMeyer(nu, g float64) float64 {
	if nu <= 0 {
		return 1
	}
	lo, hi := 1.0, 2.0
	for PrandtlMeyer(hi, g) < nu {
		lo, hi = hi, 2*hi
	}
	M := 0.5 * (lo + hi)
	for i := 0; i < 100; i++ {
		f := PrandtlMeyer(M, g) - nu
		if f > 0 {
			hi = M
		} else {
			lo = M
		}
		// dnu/dM = sqrt(M^2-1) / (M (1 + (g-1)/2 M^2))
		d := math.Sqrt(M*M-1) / (M * (1 + 0.5*(g-1)*M*M))
		next := M - f/d
		if d == 0 || next <= lo || next >= hi {
			next = 0.5 * (lo + hi)
		}
		if math.Abs(next-M) < 1e-14*M {
			return next
		}
		M = next
	}
	return M
}

// AreaRatio is the isentropic A/A* (h/h* for planar flow).
func AreaRatio(M, g float64) float64 {
	t := 2 / (g + 1) * (1 + 0.5*(g-1)*M*M)
	return math.Pow(t, (g+1)/(2*(g-1))) / M
}

// PressureRatio is the isentropic static-to-stagnation pressure p/p0.
func PressureRatio(M, g float64) float64 {
	return math.Pow(1+0.5*(g-1)*M*M, -g/(g-1))
}

// Point is a characteristic-net node.
type Point struct {
	X, Y      float64
	Theta, Nu float64 // flow angle and Prandtl-Meyer angle, rad
	M, Mu     float64
}

func newPoint(theta, nu, g float64) Point {
	M := MachFromPrandtlMeyer(nu, g)
	return Point{Theta: theta, Nu: nu, M: M, Mu: math.Asin(1 / M)}
}

// Design holds the characteristic net and wall contour.
type Design struct {
	Gamma, Me, ThetaMax float64
	N                   int
	Net                 [][]Point // Net[i][j], j >= i: C+ line i crossing C- line j
	Wall                []Point   // Wall[0] is the throat corner (state just past the fan)
}

// intersect returns the crossing of the line through a at angle ta with the
// line through b at angle tb.
func intersect(a Point, ta float64, b Point, tb float64) (float64, float64) {
	ma, mb := math.Tan(ta), math.Tan(tb)
	x := (b.Y - a.Y + ma*a.X - mb*b.X) / (ma - mb)
	return x, a.Y + ma*(x-a.X)
}

// MinimumLengthNozzle builds the nozzle with n characteristics in the corner fan.
func MinimumLengthNozzle(Me, g float64, n int) (*Design, error) {
	if Me <= 1 || g <= 1 || n < 2 {
		return nil, errors.New("need Me > 1, gamma > 1, n >= 2")
	}
	thMax := PrandtlMeyer(Me, g) / 2

	// Fan line j leaves the corner with theta_j = nu_j (K+ = 0 from the sonic
	// throat). Equal spacing; the first wave is a full step off sonic, which is
	// the usual source of the small O(1/n) area error.
	fan := make([]float64, n)
	for j := range fan {
		fan[j] = thMax * float64(j+1) / float64(n)
	}
	corner := func(j int) Point {
		p := newPoint(fan[j], fan[j], g)
		p.X, p.Y = 0, 1
		return p
	}

	d := &Design{Gamma: g, Me: Me, ThetaMax: thMax, N: n, Net: make([][]Point, n)}
	for i := 0; i < n; i++ {
		d.Net[i] = make([]Point, n)
		for j := i; j < n; j++ {
			// C- line j carries K- = 2 theta_j; the C+ line from centreline
			// point (i,i) carries K+ = -2 theta_i. Hence:
			p := newPoint(fan[j]-fan[i], fan[j]+fan[i], g)

			up := corner(j) // upstream neighbour on the C- line
			if i > 0 {
				up = d.Net[i-1][j]
			}
			cm := 0.5 * ((up.Theta - up.Mu) + (p.Theta - p.Mu))
			if j == i {
				// Centreline: theta = 0, y = 0.
				p.X, p.Y = up.X-up.Y/math.Tan(cm), 0
			} else {
				left := d.Net[i][j-1]
				cp := 0.5 * ((left.Theta + left.Mu) + (p.Theta + p.Mu))
				p.X, p.Y = intersect(left, cp, up, cm)
			}
			d.Net[i][j] = p
		}
	}

	// Wall: each C+ line, after crossing the last fan wave, runs through a
	// uniform region to the wall; the wall there takes that flow angle so the
	// wave is cancelled rather than reflected.
	w0 := corner(n - 1)
	d.Wall = append(d.Wall, w0)
	for i := 0; i < n; i++ {
		last := d.Net[i][n-1]
		prev := d.Wall[len(d.Wall)-1]
		w := last
		w.X, w.Y = intersect(last, last.Theta+last.Mu, prev, 0.5*(prev.Theta+w.Theta))
		d.Wall = append(d.Wall, w)
	}
	return d, nil
}

// ExitHeight is h_e/h*, which mass conservation says should equal A/A*(Me).
func (d *Design) ExitHeight() float64 { return d.Wall[len(d.Wall)-1].Y }

// Length is the throat-to-exit length in throat half-heights.
func (d *Design) Length() float64 { return d.Wall[len(d.Wall)-1].X }

// WallThrustCoefficient integrates axial momentum into the half-nozzle,
// normalised by p0 h* (i.e. F / (p0 A*) per unit depth):
//
//	throat     : p* (1 + g) h*      (uniform sonic, rho a^2 = g p)
//	wall force : integral of p dy   (pressure on the diverging wall)
//
// With no ambient pressure this is the vacuum thrust coefficient.
func (d *Design) WallThrustCoefficient() float64 {
	g := d.Gamma
	cf := PressureRatio(1, g) * (1 + g)
	for k := 1; k < len(d.Wall); k++ {
		a, b := d.Wall[k-1], d.Wall[k]
		cf += 0.5 * (PressureRatio(a.M, g) + PressureRatio(b.M, g)) * (b.Y - a.Y)
	}
	return cf
}

// MomentumBalance compares the throat-plus-wall momentum with the exit-plane
// momentum flux pe (1 + g Me^2) he and returns (in - out) / out. It depends on
// the wall contour and its pressures, so it tests the whole solution, not just
// the exit state.
func (d *Design) MomentumBalance() float64 {
	g := d.Gamma
	out := PressureRatio(d.Me, g) * (1 + g*d.Me*d.Me) * d.ExitHeight()
	return (d.WallThrustCoefficient() - out) / out
}

// ThrustCoefficient is F / (p0 A*) for uniform axial exit flow at Me with
// area ratio eps and ambient ratio pa/p0. Continuity gives
// mdot Ve / (p0 A*) = g Me^2 (pe/p0) eps.
func ThrustCoefficient(Me, g, eps, paRatio float64) float64 {
	pe := PressureRatio(Me, g)
	return g*Me*Me*pe*eps + (pe-paRatio)*eps
}

// WedgeThrustCoefficient applies the planar divergence factor to the
// momentum term only; the exit pressure still acts over the full exit area.
func WedgeThrustCoefficient(Me, g, eps, paRatio, alpha float64) float64 {
	pe := PressureRatio(Me, g)
	return WedgeDivergenceFactor(alpha)*g*Me*Me*pe*eps + (pe-paRatio)*eps
}

// IdealCf is the textbook closed form for an ideal 1-D nozzle (Sutton eq. 3-30).
func IdealCf(g, peRatio, eps, paRatio float64) float64 {
	k := 2 * g * g / (g - 1) * math.Pow(2/(g+1), (g+1)/(g-1))
	return math.Sqrt(k*(1-math.Pow(peRatio, (g-1)/g))) + (peRatio-paRatio)*eps
}

// WedgeDivergenceFactor is the axial-momentum fraction for a planar straight
// wedge of half-angle alpha (source flow): sin(alpha)/alpha.
func WedgeDivergenceFactor(alpha float64) float64 { return math.Sin(alpha) / alpha }

// WedgeLength is the divergent length of a planar wedge reaching half-height
// eps from a unit throat.
func WedgeLength(eps, alpha float64) float64 { return (eps - 1) / math.Tan(alpha) }
