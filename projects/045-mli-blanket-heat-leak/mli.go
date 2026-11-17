// Astron Built this AP
// Project 045 | Thermal analysis | Multilayer insulation layer-by-layer heat leak
// Author: Aditya (aditya-astron)
// Date: 2026-11-17 IST

// Layer-by-layer model of a multilayer insulation (MLI) blanket.
//
// Nodes 0..n+1: node 0 is the hot boundary (innermost layer, held at the wall
// temperature), nodes 1..n are the radiation shields, node n+1 is the cold
// boundary (either a fixed temperature or the outer cover radiating to space).
// Gap k sits between node k and node k+1 and carries, per unit area,
//
//	q_k = sigma (Ta^4 - Tb^4) / (1/ea + 1/eb - 1)      parallel grey plates
//	    + (hs/2) (Ta^2 - Tb^2)                        spacer conduction
//	    + Cg P (Ta^0.52 - Tb^0.52)                    free-molecular gas
//
// The spacer term is the Lockheed solid-conduction term applied per gap:
// a conductance hs * Tm with hs = Cs Nbar^2.63, Tm the gap mean temperature,
// gives hs Tm (Ta - Tb) = (hs/2)(Ta^2 - Tb^2). Summed over n+1 equal-flux gaps
// it telescopes to exactly Cs Nbar^2.63 (Th^2 - Tc^2) / (2 (n+1)), which is the
// Lockheed solid term. The gas term uses the Lockheed temperature exponent.
//
// Steady state: every gap carries the same q. With uniform emissivity every
// gap flux is F(Ta) - F(Tb) for a single potential F (a Kirchhoff transform),
// so q = (F(Th) - F(Tc)) / (n+1) in closed form. The Newton solver does not
// use that; it is the reference it is checked against.

package main

import (
	"errors"
	"math"
)

const sigma = 5.670374419e-8 // W m^-2 K^-4

// Lockheed empirical constants (Keller, Cunnington and Glassford 1974; DAM
// shields on silk net). Units: W/m^2 with T in K, Nbar in layers/cm, P in torr.
const (
	lockheedCr = 5.39e-10
	lockheedCs = 8.95e-8
	lockheedCg = 1.46e4
)

// Blanket describes the stack and its boundary conditions.
type Blanket struct {
	N        int       // number of shields between the two boundary layers
	EpsHot   []float64 // emissivity of the face of node i facing the hot side, len N+2
	EpsCold  []float64 // emissivity of the face of node i facing the cold side, len N+2
	Density  float64   // layer density Nbar, layers/cm
	Pressure float64   // residual gas pressure, torr
	NoRad    bool      // switch radiation off (verification only)

	Th float64 // hot boundary temperature, K

	// Cold boundary: fixed temperature if Radiating is false, otherwise an
	// outer cover with outer-face absorptance Alpha and emittance EpsOut that
	// absorbs Flux (W/m^2, sunlight normal to the blanket) and sees 0 K space.
	Tc        float64
	Radiating bool
	Alpha     float64
	EpsOut    float64
	Flux      float64
}

// UniformBlanket builds a stack with the same emissivity on every face.
func UniformBlanket(n int, eps, density, pressure, th, tc float64) *Blanket {
	b := &Blanket{N: n, Density: density, Pressure: pressure, Th: th, Tc: tc}
	b.EpsHot = make([]float64, n+2)
	b.EpsCold = make([]float64, n+2)
	for i := range b.EpsHot {
		b.EpsHot[i], b.EpsCold[i] = eps, eps
	}
	return b
}

func (b *Blanket) hs() float64 { return lockheedCs * math.Pow(b.Density, 2.63) }

// radFactor is sigma/(1/ea + 1/eb - 1) for gap k.
func (b *Blanket) radFactor(k int) float64 {
	if b.NoRad {
		return 0
	}
	ea, eb := b.EpsCold[k], b.EpsHot[k+1]
	return sigma / (1/ea + 1/eb - 1)
}

// GapFlux returns the heat flux across gap k and its partial derivatives with
// respect to the hot-side and cold-side node temperatures.
func (b *Blanket) GapFlux(k int, ta, tb float64) (q, dqa, dqb float64) {
	r := b.radFactor(k)
	h := b.hs()
	g := lockheedCg * b.Pressure
	q = r*(pow4(ta)-pow4(tb)) + 0.5*h*(ta*ta-tb*tb) + g*(math.Pow(ta, 0.52)-math.Pow(tb, 0.52))
	dqa = 4*r*ta*ta*ta + h*ta + 0.52*g*math.Pow(ta, -0.48)
	dqb = -(4*r*tb*tb*tb + h*tb + 0.52*g*math.Pow(tb, -0.48))
	return
}

func pow4(x float64) float64 { x2 := x * x; return x2 * x2 }

// Solution holds node temperatures 0..N+1 and the through-flux.
type Solution struct {
	T     []float64
	Q     float64 // flux leaving the wall into the blanket, W/m^2
	Iters int
}

// SolveNewton solves the steady nodal balance with Newton's method. The
// Jacobian is tridiagonal (each node couples to its two neighbours), so each
// step is one Thomas sweep. Unknowns: T1..TN, plus T(N+1) for a radiating cover.
func (b *Blanket) SolveNewton() (*Solution, error) {
	n := b.N
	m := n
	if b.Radiating {
		m = n + 1
	}
	T := make([]float64, n+2)
	T[0] = b.Th
	tEnd := b.Tc
	if b.Radiating {
		// Guess: cover in radiative equilibrium with the sun alone, floored.
		tEnd = math.Max(math.Pow(b.Alpha*b.Flux/(b.EpsOut*sigma), 0.25), 0.5*b.Th)
	}
	for i := 1; i <= n+1; i++ {
		T[i] = b.Th + (tEnd-b.Th)*float64(i)/float64(n+1)
	}

	lo, di, up, rhs := make([]float64, m), make([]float64, m), make([]float64, m), make([]float64, m)
	for it := 1; it <= 100; it++ {
		// Residual at node i (1-based): q_{i-1} - q_i (heat in minus heat out).
		qs := make([]float64, n+1)
		dA := make([]float64, n+1)
		dB := make([]float64, n+1)
		for k := 0; k <= n; k++ {
			qs[k], dA[k], dB[k] = b.GapFlux(k, T[k], T[k+1])
		}
		for j := 0; j < m; j++ {
			i := j + 1
			var res float64
			if i <= n {
				res = qs[i-1] - qs[i]
				lo[j] = dA[i-1]
				di[j] = dB[i-1] - dA[i]
				up[j] = -dB[i]
			} else { // radiating cover: q_n + alpha S - eps sigma T^4 = 0
				ti := T[i]
				res = qs[n] + b.Alpha*b.Flux - b.EpsOut*sigma*pow4(ti)
				lo[j] = dA[n]
				di[j] = dB[n] - 4*b.EpsOut*sigma*ti*ti*ti
				up[j] = 0
			}
			rhs[j] = -res
		}
		dx := thomas(lo, di, up, rhs)
		step := 0.0
		for j := 0; j < m; j++ {
			// Cap the step so an early overshoot cannot drive T negative.
			d := math.Max(math.Min(dx[j], 50), -0.5*T[j+1])
			T[j+1] += d
			step = math.Max(step, math.Abs(d))
		}
		if step < 1e-11*b.Th {
			q, _, _ := b.GapFlux(0, T[0], T[1])
			return &Solution{T: T, Q: q, Iters: it}, nil
		}
	}
	return nil, errors.New("newton did not converge")
}

// thomas solves a tridiagonal system (lo[0] and up[m-1] unused).
func thomas(lo, di, up, rhs []float64) []float64 {
	m := len(di)
	c := make([]float64, m)
	d := make([]float64, m)
	c[0] = up[0] / di[0]
	d[0] = rhs[0] / di[0]
	for i := 1; i < m; i++ {
		den := di[i] - lo[i]*c[i-1]
		c[i] = up[i] / den
		d[i] = (rhs[i] - lo[i]*d[i-1]) / den
	}
	x := make([]float64, m)
	x[m-1] = d[m-1]
	for i := m - 2; i >= 0; i-- {
		x[i] = d[i] - c[i]*x[i+1]
	}
	return x
}

// SolveShooting is an independent check: guess q, march from the wall solving
// each gap for the next temperature, and bisect q until the far boundary
// condition holds. It uses only GapFlux, never the Jacobian or the potential.
func (b *Blanket) SolveShooting() *Solution {
	march := func(q float64) []float64 {
		T := make([]float64, b.N+2)
		T[0] = b.Th
		for k := 0; k <= b.N; k++ {
			ta := T[k]
			// Gap flux is strictly decreasing in Tb; bisect on a wide bracket.
			lo, hi := 1e-6, 1e4
			for it := 0; it < 200; it++ {
				mid := 0.5 * (lo + hi)
				f, _, _ := b.GapFlux(k, ta, mid)
				if f > q {
					lo = mid
				} else {
					hi = mid
				}
			}
			T[k+1] = 0.5 * (lo + hi)
		}
		return T
	}
	// Mismatch is increasing in q for both boundary types.
	mismatch := func(q float64) float64 {
		T := march(q)
		te := T[b.N+1]
		if b.Radiating {
			return q + b.Alpha*b.Flux - b.EpsOut*sigma*pow4(te)
		}
		return b.Tc - te
	}
	qlo, qhi := -1e4, 1e4
	for it := 0; it < 200; it++ {
		mid := 0.5 * (qlo + qhi)
		if mismatch(mid) > 0 {
			qhi = mid
		} else {
			qlo = mid
		}
	}
	q := 0.5 * (qlo + qhi)
	return &Solution{T: march(q), Q: q}
}

// Potential is the Kirchhoff potential F(T) for a uniform-emissivity stack.
func (b *Blanket) Potential(t float64) float64 {
	r := b.radFactor(0)
	return r*pow4(t) + 0.5*b.hs()*t*t + lockheedCg*b.Pressure*math.Pow(t, 0.52)
}

// ClosedForm returns q and the node temperatures for a uniform stack with
// fixed boundaries: F drops by the same amount q across every gap.
func (b *Blanket) ClosedForm() *Solution {
	n1 := float64(b.N + 1)
	fh, fc := b.Potential(b.Th), b.Potential(b.Tc)
	q := (fh - fc) / n1
	T := make([]float64, b.N+2)
	for k := range T {
		target := fh - float64(k)*q
		lo, hi := 0.0, b.Th*2
		for it := 0; it < 200; it++ {
			mid := 0.5 * (lo + hi)
			if b.Potential(mid) > target {
				hi = mid
			} else {
				lo = mid
			}
		}
		T[k] = 0.5 * (lo + hi)
	}
	T[0], T[b.N+1] = b.Th, b.Tc
	return &Solution{T: T, Q: q}
}

// Lockheed returns the empirical Lockheed heat flux for n shields at the
// given boundary temperatures (radiation, solid and gas terms).
func Lockheed(n int, eps, density, pressure, th, tc float64) (qr, qs, qg float64) {
	nf := float64(n)
	qr = lockheedCr * eps * (math.Pow(th, 4.67) - math.Pow(tc, 4.67)) / nf
	qs = lockheedCs * math.Pow(density, 2.63) * (th - tc) * (th + tc) / (2 * (nf + 1))
	qg = lockheedCg * pressure * (math.Pow(th, 0.52) - math.Pow(tc, 0.52)) / nf
	return
}

// EffectiveEmittance converts a flux to the figure of merit
// eps* = q / (sigma (Th^4 - Tc^4)) used in spacecraft thermal models.
func EffectiveEmittance(q, th, tc float64) float64 {
	return q / (sigma * (pow4(th) - pow4(tc)))
}

// FixedThicknessFlux is the closed-form uniform-stack flux when n shields are
// packed into a blanket of thickness tcm (cm), so Nbar = n / tcm. n may be
// real so the optimum can be found by a 1-D search.
func FixedThicknessFlux(n, tcm, eps, pressure, th, tc float64) float64 {
	b := &Blanket{Density: n / tcm, Pressure: pressure}
	b.EpsHot = []float64{eps, eps}
	b.EpsCold = []float64{eps, eps}
	return (b.Potential(th) - b.Potential(tc)) / (n + 1)
}

// GoldenMin minimises a unimodal f on [a, c].
func GoldenMin(f func(float64) float64, a, c float64) float64 {
	gr := (math.Sqrt(5) - 1) / 2
	x1, x2 := c-gr*(c-a), a+gr*(c-a)
	f1, f2 := f(x1), f(x2)
	for c-a > 1e-10 {
		if f1 < f2 {
			c, x2, f2 = x2, x1, f1
			x1 = c - gr*(c-a)
			f1 = f(x1)
		} else {
			a, x1, f1 = x1, x2, f2
			x2 = a + gr*(c-a)
			f2 = f(x2)
		}
	}
	return 0.5 * (a + c)
}
