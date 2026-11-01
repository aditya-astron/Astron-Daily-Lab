# 004 - Solar Array MPPT Through a LEO Eclipse Cycle

**Domain:** Spacecraft power systems

## Problem

A sun-pointed array of triple-junction GaAs cells (16 in series x 12 strings,
~230 W at 28 C) flies a 500 km orbit with the longest eclipse (beta = 0). It
leaves every eclipse near -90 C and warms to about +40 C over the sunlit arc.
Vmp moves by ~14 V over that swing. How much energy does a fixed-voltage
(direct energy transfer) design leave on the table, how well do perturb and
observe and incremental conductance trackers do, and how should the tracker
step size be chosen?

All you get from a cell datasheet is Isc, Voc, Vmp, Imp and two temperature
coefficients, so the first job is to turn those into a model.

## Approach

**Cell model.** Single diode, per cell:

```
I = Iph - I0 (exp((V + I Rs)/a) - 1) - (V + I Rs)/Rsh,    a = n k T / q
```

- Parameters come from four conditions: `I(0) = Isc`, `I(Voc) = 0`,
  `I(Vmp) = Imp`, and `dP/dV = 0` at the MPP. For a trial `Rs` the first three
  are linear in `(Iph, I0, 1/Rsh)`, so they are solved as a 3x3 system and only
  the MPP condition `G/(1 + G Rs) = Imp/Vmp` is left, a 1-D Brent root in `Rs`.
  (`G = (I0/a) exp(Vd/a) + 1/Rsh` is the diode plus shunt conductance.)
- The ideality factor is an input. For this cell's 0.864 fill factor a
  physical fit (`Rsh > 0`) exists only for `n < ~3.05`; `n = 2.9` gives
  `Rs = 52.6 mOhm`, `Rsh = 2.82 kOhm`, `I0 = 1.36e-16 A`.
- The current is evaluated explicitly with Lambert W. Its argument is around
  `exp(400)` here, so it is computed as the Wright omega function,
  `omega(x) = W(exp(x))`, from the logarithm.
- Temperature: `Iph` follows the Isc coefficient and scales with irradiance;
  `I0` is re-derived at each temperature so Voc follows the datasheet
  coefficient exactly. The Pmp coefficient then falls out as -0.23 %/K.

**Environment.** Cylindrical Earth shadow with a 9 s linear penumbra ramp.
Panel temperature from a one-node balance per m^2:

```
C dT/dt = (alpha - eta) G + q_IR - (eps_f + eps_b) sigma T^4
```

with alpha 0.91, emissivities 0.85/0.80, 70 W/m^2 Earth IR, 2500 J/m^2/K. The
start temperature is iterated until the orbit is periodic.

**Trackers** (10 Hz, ideal converter, 20 mV / 10 mA sensor noise):

- P&O: step in the same direction if the last step raised power, else reverse.
- IncCond: at the MPP `dI/dV = -I/V`; the sign of `dI/dV + I/V` gives the
  direction, with a hold band of 0.5 % of `I/V`.
- In eclipse both park at 40 V and reset. The first lit sample is only used as
  a baseline: differencing it against the dark reading would be all
  irradiance change and no voltage change, which sends IncCond the wrong way
  half the time.
- DET clamps the array at 35 V, just under the hot Vmp.

Energy and efficiency are taken over the last of three orbits.

## Results

```
Orbit: 500 km, beta 0 deg, period 94.62 min, eclipse 35.75 min
Panel temperature over orbit: -92.7 C to 41.5 C
Vmp in full sun: 37.14 V to 51.53 V
Pmp in full sun: 225.7 W to 289.9 W
Energy available at MPP: 227.07 Wh per orbit

method     Wh/orbit  efficiency  99% settle
DET          213.08      93.84%       never
P&O          226.86      99.91%       6.0 s
IncCond      226.86      99.90%       7.2 s

P&O step size sweep (10 Hz updates)
 step V  efficiency  99% settle
   0.02     98.321%     170.6 s
   0.05     99.791%      48.8 s
   0.20     99.907%       6.0 s
   0.50     99.766%       4.8 s
   1.00     99.262%       never
   2.00     96.764%       never
```

- DET gives up 6.2 % of the orbit's energy, almost all of it in the first
  half of the sunlit arc while the panel is cold and the peak power is
  highest.
- With a sensible step both trackers are within 0.1 % of ideal; their
  difference is noise. Both overshoot Vmp by ~1 V during the penumbra ramp,
  because rising irradiance makes every step look like a gain. It costs
  little here since the ramp lasts seconds.
- The step size is the real design choice. Too small and the tracker spends
  minutes walking up to the cold Vmp; too large and the three-point
  oscillation around the peak never gets inside 1 %.

## Verification (`test_mppt.py`)

- Fit reproduces Isc, Imp, Voc to 1e-10 A and `dP/dV = 0` at the datasheet MPP.
- Hand calcs: fill factor 0.8643, efficiency 29.39 % at AM0.
- Lambert W current agrees with Newton on the implicit equation to 1e-12 A
  at -90, 28 and 80 C.
- Golden-section MPP agrees with a 400k-point sweep to 1e-9 relative.
- Voc and Isc temperature slopes match the datasheet; Pmp slope within
  -0.18 to -0.30 %/K; half-sun Voc drops by `a ln 2`.
- Orbit period by hand (94.62 min), shadow half-angle `asin(Re/r)`, no
  eclipse at beta 75 deg, umbra duration from the irradiance profile.
- Panel temperature under constant sun converges to the radiative
  equilibrium solved by hand.
- Both trackers converge to the true MPP on a static curve; full orbit run
  keeps every method at or below the MPP power.

## Run

```
pip install numpy scipy pytest
python mppt_sim.py
python -m pytest -q
```
