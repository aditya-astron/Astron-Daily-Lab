# 046 - S3R shunt regulator for a 50 V GEO bus

A GEO platform with a fully regulated 50 V bus needs 1.8 kW of load plus
0.7 kW of battery recharge in sunlight at end of life. The solar array is split
into sections, each one shorted to return by its own shunt switch or released
onto the bus. This project sizes the strings, the section count and the bus
capacitor of a sequential switching shunt regulator (S3R), then simulates the
switching limit cycle exactly and checks it against the closed-form theory.

## Approach

**Cell and string.** Single-diode cell model (three junctions lumped,
n = 3.3, Rs = 40 mohm) with Iph and I0 pinned to Isc and Voc in closed form,
and I(V) solved explicitly with the Lambert W function. Temperature and
end-of-life degradation act on Isc and Voc. A string sees the bus through a
blocking diode and 0.6 ohm of harness. The string length is the shortest that
keeps the bus on the flat part of the curve at EOL, 70 C:
I(50 V) >= 98 % of Isc. That gives 27 cells, and the section current then
moves by only 0.6 mA across the comparator band, so a section is a clean
current source.

**Section count.** Worst section current is at EOL, summer solstice
(cos 23.44 deg / 1.0167^2 = 0.888 sun). Required 50 A / 1.763 A rounded up,
plus one failed section: 30 sections.

**Limit cycle.** Section k shunts when V_bus reaches hi_k and releases at
lo_k = hi_k - h. Window centres are spaced s = 50 mV apart with h = 40 mV, so
exactly one section switches. With the bus capacitor C and its ESR,

```
C dVc/dt = i_c,           V_bus = Vc + ESR i_c
h_eff    = h - ESR Is     (the ESR step eats part of the band)
f        = x (Is - x) / (C h_eff Is),    f_max = Is / (4 C h_eff) at x = Is/2
duty     = x / Is         (charge balance on the active section)
mean V   = centre_j + ESR (x - Is/2)
```

where x is the load current not covered by the fully-on sections. With
constant currents every segment is a straight line, so the simulator jumps
from one comparator crossing to the next with no time step. C is sized from
f_max <= 5 kHz at the largest section current (BOL, equinox, 70 C).

**Checks** (`test_s3r.py`, 29 tests)

- Lambert W current against a bracketed root of the implicit diode equation
  to 1e-12 A, including the overflow branch far past Voc; I(0) = Isc and
  I(Voc) = 0 exactly; Vmp 2.39 V, FF 0.855, in line with commercial 3J cells.
- Event-driven frequency, duty and exact time-averaged bus voltage against the
  formulas above to 1e-9 for three ESR values and five load fractions; only
  the active section ever switches and V_bus stays inside its window.
- f_max sits at x = Is/2 and equals Is / (4 C h_eff).
- Static regulation staircase: mean bus voltage vs load has slope -s/Is.
- 20 A load step: minimum bus voltage lands exactly on the lower edge of the
  new active window, overshoot on release on the upper edge of the old one,
  and the original limit cycle returns.
- Charge balance C dVc = integral of i_c over a two-step transient.
- Independent fixed-step Euler integration with sampled comparators agrees
  on frequency to 0.2 %.
- Constructor rejects ESR Is >= h (chatter) and overlapping windows.

## Run

```
pip install numpy scipy pytest
python3 run_s3r.py
python3 -m pytest -q
```

## Sample output

```
String: 27 cells in series, 4 strings per section (knee rule 98% of Isc at EOL 70 C)
Required bus current 50.0 A, worst section 1.763 A -> 30 sections (N+1)
BOL equinox margin: 13.7 A shunted at design load (6.5 sections)
Section current change across the 40 mV band, EOL hot: 0.57 mA (3.2e-04 of Is)

Bus capacitor: h_eff = 31.50 mV -> C = 3.40 mF, f_max = 4960 Hz at Is = 2.125 A
Regulation band: 49.275 .. 50.725 V, sunlight output resistance s/Is = 23.5 mohm

   I_load   win    x/Is   f_sim Hz  f_theory Hz   duty   Vmean sim   Vmean th
    32.04    14    0.08     1460.1       1460.1   0.080    49.97143   49.97143
    32.93    14    0.50     4959.6       4959.6   0.500    49.97500   49.97500
    33.83    14    0.92     1460.1       1460.1   0.920    49.97857   49.97857

20 A load step from 22.3 A: active window 19 -> 10
  min V_bus 49.7550 V (lo of window 10: 49.7550 V)
  max V_bus after release 50.2450 V (hi of window 19: 50.2450 V)
```

## Limits

Window thresholds are bus-referred; in hardware they sit on the main error
amplifier output and are divided by its gain. Section currents are held
constant within a run (justified by the 3e-4 slope above), and the array's
own capacitance, which is dumped into the shunt switch on every shunt event,
is not modelled.
