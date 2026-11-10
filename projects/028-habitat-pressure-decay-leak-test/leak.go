// Astron Built this AP
// Project 028 | Analog habitat systems | Habitat pressure-decay leak test
// Author: Aditya (aditya-astron)
// Date: 2026-11-10 IST
//
// Pressure-decay leak test of a sealed analog habitat module. The module is
// pumped to a few kPa above ambient, sealed, and logged for several days.
// The difficulty is that the leak is small compared with what the weather
// does to the pressure: at 101 kPa and 293 K one kelvin of air temperature is
// worth p/T = 346 Pa, and the barometer moves by 100-400 Pa a day. Reading
// the leak off the raw differential pressure is therefore badly biased.
//
// The fix is to track the gas inventory rather than the pressure,
//
//	m = p_i V / (R T)                                  (ideal gas, rigid hull)
//
// and fit a leak law to it with the measured temperature and barometer as
// inputs. The leak is the building-envelope power law (ASTM E779 form),
//
//	Q = C |dp|^n sign(dp),   dp = p_i - p_e,   0.5 <= n <= 1
//
// n = 0.5 is a sharp orifice, n = 1 laminar flow through cracks and seals;
// real hulls land in between. Mass balance with the volumetric flow taken at
// the upstream density:
//
//	dm/dt = -rho_up Q,  rho_up = m/V (outflow) or p_e/(R T_e) (inflow)
//
// The fit is output-error: integrate that ODE through the measured inputs and
// adjust (ln C, n, m0) by Levenberg-Marquardt so the model inventory matches
// the measured one.
package habitatleak

import (
	"errors"
	"math"
	"sort"
)

// RAir is the specific gas constant of dry air [J/(kg K)].
const RAir = 287.05

// Habitat holds the fixed properties of the module under test.
type Habitat struct {
	Volume float64 // free gas volume [m^3]
	TExt   float64 // exterior air temperature, used only for inflow density [K]
}

// LeakLaw is the power-law leak Q = C |dp|^n sign(dp) [m^3/s].
type LeakLaw struct {
	C float64 // flow coefficient [m^3/(s Pa^n)]
	N float64 // flow exponent [-]
}

// Flow returns the volumetric leak flow for a pressure difference dp [Pa].
func (l LeakLaw) Flow(dp float64) float64 {
	if dp == 0 {
		return 0
	}
	return math.Copysign(l.C*math.Pow(math.Abs(dp), l.N), dp)
}

// Series is a sampled log: interior differential pressure, barometer, and
// volume-averaged air temperature at times T.
type Series struct {
	T  []float64 // time [s], strictly increasing
	Dp []float64 // p_i - p_e [Pa]
	Pe []float64 // exterior (barometric) pressure [Pa]
	Ta []float64 // interior mean air temperature [K]
}

// Mass converts the log to gas inventory, m = (p_e + dp) V / (R T).
func (h Habitat) Mass(s Series) []float64 {
	m := make([]float64, len(s.T))
	for k := range s.T {
		m[k] = (s.Pe[k] + s.Dp[k]) * h.Volume / (RAir * s.Ta[k])
	}
	return m
}

// Inputs supplies the exterior pressure and interior temperature at any time.
type Inputs interface {
	At(t float64) (pe, ta float64)
}

// MassRate is dm/dt for inventory m at time t.
func (h Habitat) MassRate(law LeakLaw, m, pe, ta float64) float64 {
	pi := m * RAir * ta / h.Volume
	q := law.Flow(pi - pe)
	if q >= 0 {
		return -(m / h.Volume) * q
	}
	return -(pe / (RAir * h.TExt)) * q
}

// Integrate advances the inventory with classical RK4 from m0 at times[0],
// taking sub steps per interval, and returns m at every entry of times.
func (h Habitat) Integrate(law LeakLaw, in Inputs, m0 float64, times []float64, sub int) []float64 {
	out := make([]float64, len(times))
	out[0] = m0
	m := m0
	f := func(t, m float64) float64 {
		pe, ta := in.At(t)
		return h.MassRate(law, m, pe, ta)
	}
	for k := 1; k < len(times); k++ {
		dt := (times[k] - times[k-1]) / float64(sub)
		t := times[k-1]
		for j := 0; j < sub; j++ {
			k1 := f(t, m)
			k2 := f(t+dt/2, m+dt/2*k1)
			k3 := f(t+dt/2, m+dt/2*k2)
			k4 := f(t+dt, m+dt*k3)
			m += dt / 6 * (k1 + 2*k2 + 2*k3 + k4)
			t += dt
		}
		out[k] = m
	}
	return out
}

// Sampled linearly interpolates logged inputs, which is what the fit has.
type Sampled struct {
	T, Pe, Ta []float64
}

// At implements Inputs by linear interpolation (clamped at the ends).
func (s Sampled) At(t float64) (float64, float64) {
	n := len(s.T)
	if t <= s.T[0] {
		return s.Pe[0], s.Ta[0]
	}
	if t >= s.T[n-1] {
		return s.Pe[n-1], s.Ta[n-1]
	}
	k := sort.SearchFloat64s(s.T, t) // first index with T[k] >= t
	a := (t - s.T[k-1]) / (s.T[k] - s.T[k-1])
	return s.Pe[k-1] + a*(s.Pe[k]-s.Pe[k-1]), s.Ta[k-1] + a*(s.Ta[k]-s.Ta[k-1])
}

// FitResult is the output of FitLeak.
type FitResult struct {
	Law   LeakLaw
	M0    float64
	Cov   [3][3]float64 // covariance of (ln C, n, m0)
	Sigma float64       // residual standard deviation [kg]
	Iter  int
}

// FitLeak estimates (C, n, m0) from a log by Levenberg-Marquardt on the
// inventory residuals. Fitting ln C keeps C positive and makes the problem
// far better scaled than fitting C itself. guess.N is clamped to [0.5, 1].
func (h Habitat) FitLeak(s Series, guess LeakLaw, sub int) (FitResult, error) {
	meas := h.Mass(s)
	in := Sampled{T: s.T, Pe: s.Pe, Ta: s.Ta}
	model := func(p [3]float64) []float64 {
		return h.Integrate(LeakLaw{C: math.Exp(p[0]), N: p[1]}, in, p[2], s.T, sub)
	}
	resid := func(p [3]float64) ([]float64, float64) {
		mm := model(p)
		r := make([]float64, len(mm))
		ss := 0.0
		for k := range mm {
			r[k] = mm[k] - meas[k]
			ss += r[k] * r[k]
		}
		return r, ss
	}
	jac := func(p [3]float64) [][3]float64 {
		steps := [3]float64{1e-6, 1e-6, 1e-6 * p[2]}
		J := make([][3]float64, len(s.T))
		for i := 0; i < 3; i++ {
			pp, pm := p, p
			pp[i] += steps[i]
			pm[i] -= steps[i]
			fp, fm := model(pp), model(pm)
			for k := range J {
				J[k][i] = (fp[k] - fm[k]) / (2 * steps[i])
			}
		}
		return J
	}

	p := [3]float64{math.Log(guess.C), math.Min(math.Max(guess.N, 0.5), 1), meas[0]}
	r, ss := resid(p)
	lambda := 1e-3
	iters := 0
	var J [][3]float64
	for iters < 200 {
		iters++
		J = jac(p)
		var A [3][3]float64
		var g [3]float64
		for k := range J {
			for i := 0; i < 3; i++ {
				g[i] += J[k][i] * r[k]
				for j := 0; j < 3; j++ {
					A[i][j] += J[k][i] * J[k][j]
				}
			}
		}
		improved, converged := false, false
		for tries := 0; tries < 30; tries++ {
			var B [3][3]float64
			for i := 0; i < 3; i++ {
				B[i] = A[i]
				B[i][i] *= 1 + lambda
			}
			var negG [3]float64
			for i := range g {
				negG[i] = -g[i]
			}
			d, err := solve3(B, negG)
			if err != nil {
				lambda *= 10
				continue
			}
			q := [3]float64{p[0] + d[0], p[1] + d[1], p[2] + d[2]}
			rq, ssq := resid(q)
			if ssq < ss {
				rel := (ss - ssq) / ss
				p, r, ss = q, rq, ssq
				lambda = math.Max(lambda/10, 1e-12)
				improved = true
				converged = rel < 1e-12
				break
			}
			lambda *= 10
		}
		// No downhill step at any damping means we sit at the minimum to
		// within the finite-difference Jacobian; that is also convergence.
		if !improved || converged {
			break
		}
	}

	J = jac(p)
	var A [3][3]float64
	for k := range J {
		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				A[i][j] += J[k][i] * J[k][j]
			}
		}
	}
	inv, err := inv3(A)
	if err != nil {
		return FitResult{}, err
	}
	dof := float64(len(s.T) - 3)
	s2 := ss / dof
	var cov [3][3]float64
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			cov[i][j] = s2 * inv[i][j]
		}
	}
	return FitResult{
		Law:   LeakLaw{C: math.Exp(p[0]), N: p[1]},
		M0:    p[2],
		Cov:   cov,
		Sigma: math.Sqrt(s2),
		Iter:  iters,
	}, nil
}

// FlowAt returns Q(dpRef) and its 1-sigma from the fit covariance by the
// delta method: ln Q = ln C + n ln dpRef is linear in the fitted parameters.
func (f FitResult) FlowAt(dpRef float64) (q, sigma float64) {
	q = f.Law.Flow(dpRef)
	l := math.Log(dpRef)
	v := f.Cov[0][0] + l*l*f.Cov[1][1] + 2*l*f.Cov[0][1]
	return q, q * math.Sqrt(v)
}

// EquivalentArea is the sharp-edged orifice area that passes flow q at dp:
// q = Cd A sqrt(2 dp / rho). Cd = 0.61 by convention.
func EquivalentArea(q, dp, rho, cd float64) float64 {
	return q / (cd * math.Sqrt(2*dp/rho))
}

// NaiveDecay is what a test without temperature or barometric correction
// reports: fit dp(t) to the isothermal, constant-ambient decay of a pure
// power law with the module density frozen, which for n != 1 integrates to
//
//	dp^(1-n) = dp0^(1-n) - (1-n) k t,   k = p_i0 C / V
//
// With n fixed (the usual practice when only the pressure trace is used) the
// right-hand side is a straight line in t, so k comes from ordinary least
// squares on dp^(1-n). Only samples with dp > 0 are used.
func NaiveDecay(s Series, h Habitat, n float64) LeakLaw {
	var sx, sy, sxx, sxy, cnt float64
	for k := range s.T {
		if s.Dp[k] <= 0 {
			continue
		}
		x, y := s.T[k], math.Pow(s.Dp[k], 1-n)
		sx += x
		sy += y
		sxx += x * x
		sxy += x * y
		cnt++
	}
	slope := (cnt*sxy - sx*sy) / (cnt*sxx - sx*sx)
	k := -slope / (1 - n)
	pi0 := s.Pe[0] + s.Dp[0]
	return LeakLaw{C: k * h.Volume / pi0, N: n}
}

// EndpointLeak is the textbook temperature-compensated two-point method:
// mean mass loss rate between samples a and b, attributed to the mean
// differential pressure over the window.
func EndpointLeak(s Series, h Habitat, a, b int) (mdot, dpMean float64) {
	m := h.Mass(s)
	mdot = (m[a] - m[b]) / (s.T[b] - s.T[a])
	for k := a; k <= b; k++ {
		dpMean += s.Dp[k]
	}
	return mdot, dpMean / float64(b-a+1)
}

func solve3(A [3][3]float64, b [3]float64) ([3]float64, error) {
	inv, err := inv3(A)
	if err != nil {
		return [3]float64{}, err
	}
	var x [3]float64
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			x[i] += inv[i][j] * b[j]
		}
	}
	return x, nil
}

// inv3 inverts a 3x3 matrix by cofactors; the normal matrices here are
// symmetric positive definite and only 3x3, so this is exact enough.
func inv3(a [3][3]float64) ([3][3]float64, error) {
	var c [3][3]float64
	c[0][0] = a[1][1]*a[2][2] - a[1][2]*a[2][1]
	c[0][1] = a[0][2]*a[2][1] - a[0][1]*a[2][2]
	c[0][2] = a[0][1]*a[1][2] - a[0][2]*a[1][1]
	c[1][0] = a[1][2]*a[2][0] - a[1][0]*a[2][2]
	c[1][1] = a[0][0]*a[2][2] - a[0][2]*a[2][0]
	c[1][2] = a[0][2]*a[1][0] - a[0][0]*a[1][2]
	c[2][0] = a[1][0]*a[2][1] - a[1][1]*a[2][0]
	c[2][1] = a[0][1]*a[2][0] - a[0][0]*a[2][1]
	c[2][2] = a[0][0]*a[1][1] - a[0][1]*a[1][0]
	det := a[0][0]*c[0][0] + a[0][1]*c[1][0] + a[0][2]*c[2][0]
	scale := math.Abs(a[0][0]*a[1][1]*a[2][2]) + 1e-300
	if math.Abs(det) < 1e-14*scale {
		return c, errors.New("singular normal matrix")
	}
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			c[i][j] /= det
		}
	}
	return c, nil
}
