// Astron Built this AP
// Project 028 | Analog habitat systems | Habitat pressure-decay leak test
// Author: Aditya (aditya-astron)
// Date: 2026-11-10 IST
//
// Reference test scenario: a 300 m^3 rigid analog habitat module sealed at
// +3 kPa and logged for 72 h every 5 min. The truth is integrated finely
// through smooth weather; the log adds instrument noise on top.
package habitatleak

import (
	"math"
	"math/rand"
)

const day = 86400.0

// Weather is the true interior air temperature and barometric pressure.
// Interior: 293.15 K with a 1.0 K diurnal swing (HVAC deadband, solar load)
// and a 0.3 K second harmonic. Exterior: the 12 h atmospheric tide
// (~80 Pa) plus a 400 Pa pressure drop as a front passes on day 2.
type Weather struct{}

// At implements Inputs.
func (Weather) At(t float64) (pe, ta float64) {
	w := 2 * math.Pi / day
	ta = 293.15 + 1.0*math.Sin(w*t) + 0.3*math.Sin(2*w*t+1.0)
	pe = 101325 + 80*math.Sin(2*w*t+0.4) - 200*(1+math.Tanh((t-1.5*day)/(4*3600)))
	return pe, ta
}

// Instruments are 1-sigma white-noise levels of the logged channels.
type Instruments struct {
	Dp float64 // differential transducer [Pa]
	Pe float64 // barometer [Pa]
	Ta float64 // average of the air temperature sensor array [K]
}

// Scenario bundles everything that defines a test.
type Scenario struct {
	Hab      Habitat
	Truth    LeakLaw
	Dp0      float64 // sealed-off differential pressure [Pa]
	Duration float64 // [s]
	Step     float64 // log interval [s]
	Inst     Instruments
}

// Reference returns the scenario used by the demo and the tests.
func Reference() Scenario {
	return Scenario{
		Hab:      Habitat{Volume: 300, TExt: 288.15},
		Truth:    LeakLaw{C: 2.3e-7, N: 0.62},
		Dp0:      3000,
		Duration: 3 * day,
		Step:     300,
		Inst:     Instruments{Dp: 0.5, Pe: 2.0, Ta: 0.02},
	}
}

// Times returns the logging instants.
func (sc Scenario) Times() []float64 {
	n := int(math.Round(sc.Duration/sc.Step)) + 1
	t := make([]float64, n)
	for k := range t {
		t[k] = float64(k) * sc.Step
	}
	return t
}

// TrueLog integrates the truth with 10 s RK4 steps (30 per log interval)
// and returns the noiseless log plus the true inventory.
func (sc Scenario) TrueLog() (Series, []float64) {
	t := sc.Times()
	var w Weather
	pe0, ta0 := w.At(0)
	m0 := (pe0 + sc.Dp0) * sc.Hab.Volume / (RAir * ta0)
	m := sc.Hab.Integrate(sc.Truth, w, m0, t, int(sc.Step/10))
	s := Series{T: t, Dp: make([]float64, len(t)), Pe: make([]float64, len(t)), Ta: make([]float64, len(t))}
	for k := range t {
		pe, ta := w.At(t[k])
		s.Pe[k], s.Ta[k] = pe, ta
		s.Dp[k] = m[k]*RAir*ta/sc.Hab.Volume - pe
	}
	return s, m
}

// NoisyLog adds instrument noise to a clean log.
func (sc Scenario) NoisyLog(clean Series, rng *rand.Rand) Series {
	n := len(clean.T)
	s := Series{T: clean.T, Dp: make([]float64, n), Pe: make([]float64, n), Ta: make([]float64, n)}
	for k := 0; k < n; k++ {
		s.Dp[k] = clean.Dp[k] + sc.Inst.Dp*rng.NormFloat64()
		s.Pe[k] = clean.Pe[k] + sc.Inst.Pe*rng.NormFloat64()
		s.Ta[k] = clean.Ta[k] + sc.Inst.Ta*rng.NormFloat64()
	}
	return s
}
