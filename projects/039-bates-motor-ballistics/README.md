# 039 - BATES Solid Motor Internal Ballistics: Throat Sizing, Fill, Burn and Tail-off

**Domain:** Rocket propulsion

## Problem

A 98 mm motor case is loaded with four BATES segments (cylindrical grains
with a bonded outer wall and burning core and end faces) of a generic
AP/HTPB composite. The case is rated for a 5 MPa maximum expected operating
pressure. We want the throat diameter that puts peak chamber pressure
exactly at that limit, the full pressure and thrust history from igniter
pressure through tail-off, and the delivered impulse. We also want enough
independent checks that the numbers can be trusted.

| item | value |
|---|---|
| segments | 4 x D 98 mm, core 35 mm, L 164.5 mm |
| burn rate | r = a Pc^n, 7 mm/s at 6.9 MPa, n = 0.35 |
| products | gamma 1.21, M 24 g/mol, Tc 3000 K, c* efficiency 0.95 |
| density | 1750 kg/m^3 |

## Approach

**Grain geometry.** With web burned x, core diameter d+2x and segment
length L-2x:

    Ab(x) = N [ pi (d+2x)(L-2x) + (pi/2)(D^2 - (d+2x)^2) ]

Ab is a concave quadratic in x with a peak at x* = (L-2d)/6. Setting
Ab(0) = Ab(web) gives the usual "neutral" segment length L = (3D+d)/2. A
brute-force scan over L confirms that this length minimises the max/min area
ratio. Even so, a BATES grain is still +15 % progressive at mid-web.

**Chamber mass balance.** Gas mass is m = Pc V / (R Tc) and the free volume
grows as V' = Ab r, so

    dPc/dt = (R Tc / V) [ Ab r (rho_p - rho_g) - Pc At / c* ]
    dx/dt  = r = a Pc^n

The rho_g term is the gas that fills the cavity left by the regressing
surface. It lowers equilibrium pressure by (1 - rho_g/rho_p)^(1/(1-n)),
which is 0.42 % here. Setting the derivative to zero gives the equilibrium:

    Pc = (a c* Kn (rho_p - rho_g))^(1/(1-n)),   Kn = Ab / At

**Throat sizing.** rho_g depends only on the target pressure, so the throat
follows in closed form from the peak area: At = Ab_max / Kn_max.

**Nozzle.** Ideal isentropic flow at fixed gamma. The exit Mach number is
found by bisection on the supersonic branch of A/A*, and CF is the momentum
term plus (pe - pa) eps. The expansion ratio is chosen so that pe = pa at
95 % of MEOP.

**Integration.** Fixed-step RK4 on (x, Pc, m_nozzle, I, I_vac) with a
10 us step. Web burnout is located inside the step by bisection on the
step length, so the jump in burn area to zero is never straddled. The run
stops once Pc falls below 2 pa, where the nozzle unchokes.

## Verification (`make test`, 27 checks)

- Isentropic: A/A* = 1.6875 and p/p0 = 0.1278 at M = 2, gamma 1.4 (exact and
  NACA 1135); M recovered from A/A*; CF is maximised at pe = pa; vacuum CF
  rises monotonically toward its infinite-expansion limit.
- c* = 1567.4 m/s ideal, matching the hand calculation.
- Geometry: the integral of Ab dx over the web equals the loaded propellant
  volume to 1e-15; Ab = -dV/dx; peak at (L-2d)/6; the brute-force optimum L
  equals (3D+d)/2.
- Equilibrium at the sized throat is 5 MPa to 1e-15. The simulated peak is
  5 MPa to 6e-6.
- Mid-burn the raw quasi-steady error reaches 0.4 % near burnout. That is
  real fill lag in a 5 L chamber, not error. Predicting it with
  Pc = Pc_eq - tau dPc_eq/dt, where tau is the linearised fill time
  constant, brings the residue down to 4e-5.
- Ignition: a small offset from equilibrium decays to e^-1 after the
  analytic time constant tau = V / (R Tc ((1-n) At/c* + Ab r/(R Tc))),
  within 0.13 %.
- Tail-off: after burnout Pc follows exp(-t R Tc At / (V c*)) to 1e-10.
- Conservation: throat mass plus residual chamber gas equals the loaded
  propellant to 1e-14. Vacuum impulse equals CF_vac c* m_nozzle exactly.
- Step convergence: impulse and burnout time change by less than 1e-10
  between 10 us and 2.5 us steps.

## Run

    make            # builds ./motor and ./test_ballistics (C11, libm)
    ./motor         # design summary
    ./motor --csv   # t, web, Pc, thrust, Ab trace every 20 ms
    make test

## Sample output

    BATES motor, 4 x (D 98.0 mm, d 35.0 mm, L 164.5 mm)
      web 31.50 mm, propellant 7.578 kg, volumetric loading 0.859
      Ab start/peak/end  1250.0 / 1437.0 / 1250.0 cm^2  (peak at x = 15.75 mm, +14.96 %)
    Gas: c* 1488.9 m/s (ideal 1567.3), r(5 MPa) 6.25 mm/s
    Nozzle: throat 24.39 mm dia, Kn 268..308, eps 6.60, pe/pc 0.0213
      CF sea level at MEOP 1.5482, CF vacuum 1.6820
    Simulation (dt 10 us, ignition at 0.3 MPa)
      peak Pc 5.0000 MPa (target 5.0), peak thrust 3615 N
      web burnout 5.162 s at 4.054 MPa, run ends 5.209 s (Pc < 2 pa)
      total impulse  sea level 17349 N s, vacuum 18975 N s
      delivered Isp  sea level 233.5 s, vacuum 255.3 s
      mass balance   nozzle 7.57707 + chamber 0.00078 = 7.57785 kg vs 7.57785 kg loaded
      impulse class  N (avg thrust 3331 N)

## Limits

This is a lumped (0-D) chamber model. It leaves out erosive burning (worth
checking here, since the port-to-throat area ratio starts near 2), throat
erosion, nozzle flow separation during start and tail-off (thrust is simply
clipped at zero), two-phase and frozen-flow losses beyond the single c*
efficiency, and burn-rate temperature sensitivity.
