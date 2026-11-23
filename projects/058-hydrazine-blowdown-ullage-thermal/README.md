# 058 - Hydrazine blowdown system sizing with ullage-gas thermodynamics

## Problem

A 180 kg smallsat needs 120 m/s from four 1 N hydrazine thrusters fed by a
single spherical diaphragm tank in blowdown (no regulator). The thrusters are
qualified for 5.5 to 24 bar feed. Size the tank: blowdown ratio, ullage,
loaded pressure and shell. Then answer the question first-pass sizing skips:
the N2 ullage cools as it expands, so does the feed pressure still hold
5.5 bar at depletion, and how much extra pressurant does it take?

## Approach

**Thruster.** Feed, valve, injector and catalyst-bed losses are lumped into
one orifice, `mdot = K sqrt(Pf - Pc)`, with a choked throat
`mdot = Pc At / c*`. Equating them gives a quadratic in Pc, solved in its
cancellation-free form. `At` and `K` come from the nominal point
(1 N, 22 bar feed, 11 bar chamber). Isp = c* Cf / g0 = 223.4 s, held constant.

**Blowdown sizing.** With loaded volume Vl, usable volume Vu and blowdown
ratio B = P0/Pend under `P V^n = const`:

    Vg0 = Vu / (B^(1/n) - 1),    Vt = Vl + Vg0

A pressure-sized sphere has mass proportional to `P0 Vt`. Minimising it gives
two closed forms:

    n = 1:              B* = 1 + sqrt(Vu/Vl)   (2 with no residual)
    Vl = Vu, any n:     B* = ((n+1)/n)^n       (2.127 for n = 1.4, tends to e)

The shell is Ti-6Al-4V, 1.5 burst factor, with a 0.6 mm minimum gauge.

**Ullage energy balance.** The gas is a closed system whose volume grows as
propellant leaves through the diaphragm:

    mg cv dT/dt = -P dVg/dt + h A_g (Tw - T),    P = mg R T / Vg

The wall and propellant are held at Tw: they carry about 40 kJ/K against
about 60 J/K for the gas. Wetted area scales with ullage fraction. RK4 runs the
manoeuvre plan (900 s burns, 5400 s coasts), and the depletion instant is
found inside the step by bisection. A bisection on P0 (or on Vg0 when P0 hits
the 24 bar cap) restores the floor.

## Findings

- The classic `P0 Vt` optimum (B = 1.985 here) does not apply. At about 13 to
  29 L and 8 to 24 bar, the membrane thickness is 0.13 to 0.29 mm, so every
  candidate is gauge-limited. The shell then scales with Vt^(2/3), so the
  smallest tank wins: the highest B the thrusters accept. B = 4 gives
  22 bar, 13.0 L and 0.71 kg, against 0.93 kg at B = 2.
- The thermal penalty depends on how long expulsion takes compared with the
  ullage time constant τ = mg cv / (h A_g), which is 183 s at h = 5 W/m²K. It
  does not depend on the burn plan. With 1 N thrusters, expulsion takes about
  10,000 s, so the gas trails the wall by only `P Qv / (h A_g)`, about 0.3 K.
  It dips to 285.6 K in the first burn and the floor falls just 4 mbar.
  Holding 5.5 bar takes +0.02 bar of loaded pressure, or +0.05 bar at h = 1.
- Things change when the same tank feeds bigger thrusters. With 50 N per
  thruster, expulsion takes 221 s, which is about τ. The floor falls to
  5.13 bar, the gas reaches 238 K, and P0 has to rise to 23.9 bar, right at
  the thruster limit.
- The fully adiabatic bound (h = 0) would leave 3.16 bar at depletion. No P0
  under 24 bar fixes that, so treating the ullage as adiabatic during design
  would push B well below 4 for no real reason.

## Verification (`make test`, 29 checks)

- The thruster passes through its nominal point. The Pc root matches a
  bisection on the flow balance to 1e-16, and both asymptotes are checked
  (Pc tends to Pf at low feed and to `K c* sqrt(Pf)/At` at high feed).
- The propellant budget matches a hand calculation (9.5925 kg). Both
  closed-form B* expressions match golden-section search to 1e-7 for
  n = 1, 1.2, 1.4 and 1.67.
- Tank radius and membrane thickness match hand numbers (146.0 mm,
  0.269 mm).
- **Isothermal case:** the ODE depletion time matches adaptive Simpson
  quadrature of `rho dVg / (N mdot)` to 1e-15, and the end pressure is P0/B.
  Expelled mass equals rho dVg, and the integrated dV equals both
  `Isp g0 ln(m0/mf)` and the 120 m/s requirement.
- **Adiabatic case:** `P Vg^gamma` and `T Vg^(gamma-1)` are invariant to
  1e-14 across all 15 burns and coasts. The floor is 22 bar / 4^1.4.
- **Coast:** temperature recovery matches `Tw + (T0-Tw) e^(-t/tau)` to 1e-13.
- **Quasi-steady lag:** at depletion the lag matches `P Qv / (h A_g)` to
  0.7 %.
- **Large h:** h = 2000 reproduces the isothermal firing time to 1e-5.
- **Re-sizing:** the re-sized P0 and the re-sized ullage both hold exactly
  5.5 bar.

## Run

    make test     # verification
    make run      # full study, output as in sample_output.txt

Requires a C++17 compiler; no other dependencies. Takes about 10 s in total.

## Sample output (excerpt)

    == Loaded pressure that holds the 5.5 bar floor (same tank) ==
       h W/m2K     P0 bar    dP0 bar       eff. B
          50.0     22.002      0.002        4.000
           5.0     22.017      0.017        4.003
           1.0     22.048      0.048        4.009
           0.0         --  exceeds 24 bar thruster limit

     F/thr N    t_exp s   Pmin bar     Tmin K     P0 bar      Vg0 L Vt L @22bar
           1       9988      5.495      285.6     22.021      3.189       13.04
          20        526      5.381      253.3     22.511      3.283       13.13
          50        221      5.133      238.3     23.898      3.547       13.40

## Assumptions

The thruster is generic, not a specific product. Isp is held constant through
blowdown. N2 is treated as calorically perfect. The tank wall stays at a fixed
temperature, and feed-line pressure drop is ignored. The convective
coefficient h is a parameter, because microgravity ullage heat transfer is
mostly conduction and depends on the tank.
