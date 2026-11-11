# 032 - Edelbaum Low-Thrust LEO-to-GEO Transfer with Plane Change

**Domain:** Orbital mechanics

## Problem

An all-electric GEO satellite is dropped by the launcher into a 200 km,
28.5 deg circular parking orbit and has to raise itself to GEO and remove
the full 28.5 deg of inclination with a 2 N, 1800 s Isp thruster. How much
delta-V and propellant does that take, how should the thrust be yawed out of
plane along the way, and how much worse is it than an impulsive transfer or
than doing the spiral and the plane change one after the other? The steering
law also has to be shown to actually deliver the spacecraft to GEO when flown
in full two-body dynamics, not just in the averaged model it came from.

## Approach

**Averaged dynamics.** For a quasi-circular orbit of speed `v` with thrust
acceleration `f` yawed by `beta` out of the velocity direction, and the
normal component flipped at the antinodes (sign of `cos u`):

    dv/dt = -f cos(beta)
    di/dt = -(2/pi) f sin(beta) / v

**Edelbaum solution.** Minimising delta-V gives `v sin(beta) = const`. With
`V` the accumulated delta-V:

    v(V)    = sqrt(v0^2 - 2 v0 V cos(beta0) + V^2)
    beta(V) = atan2(v0 sin(beta0), v0 cos(beta0) - V)
    di(V)   = (2/pi) (beta(V) - beta0)

    tan(beta0) = sin(pi di/2) / (v0/v1 - cos(pi di/2))
    dV         = sqrt(v0^2 - 2 v0 v1 cos(pi di/2) + v1^2)

The yaw is scheduled on accumulated delta-V, so the same law stays optimal
with constant thrust and falling mass; only the clock changes, and flight
time follows from the rocket equation.

**Checks.**

- Limits: `di = 0` gives `v0 - v1`; equal radii give `2 v0 sin(pi di/4)`,
  which beats a pure out-of-plane burn `(pi/2) v0 di`.
- The averaged ODEs integrated under the yaw law (DOP853, rtol 1e-12)
  land on `v1` and `di` to 1e-9.
- Direct optimisation: 30 to 40 constant-yaw arcs (exact on each arc, since
  `v` is linear in `V`) solved by SLSQP from a deliberately wrong start never
  beats the closed form and converges to within 0.05 %, including a lowering
  case where `beta > 90 deg`.
- Baselines: sequential spiral-then-plane-change, and an impulsive Hohmann
  with the plane change optimally split between perigee and apogee.
- Cartesian flight: the steered spacecraft is propagated in two-body
  dynamics with constant thrust and mass flow, stopping when the accumulated
  delta-V reaches the Edelbaum value. Final `a`, `i` and `e` are compared
  with the target.

**Limits.** Two-body only (no J2, drag, eclipses or third bodies), continuous
thrust, quasi-circular orbits. The residual eccentricity is the forced
response to in-plane thrust, of order `f / (n^2 r)`, and shrinks with thrust.

## How to run

    pip install numpy scipy pytest
    python3 edelbaum.py
    python3 -m pytest -q

The demo takes about 25 s (the 59-day Cartesian flight), the tests about
25 s.

## Sample output

    Edelbaum LEO 200 km / 28.5 deg -> GEO / 0 deg
      v0 = 7.7843 km/s, v1 = 3.0747 km/s
      dV = 6.0053 km/s, beta0 = 21.135 deg, beta_end = 65.902 deg

      Yaw and inclination schedule (accumulated dV):
         V [km/s]  r [km]    beta [deg]  i [deg]
           0.000      6578     21.135    28.500
           1.201      8938     24.853    26.133
           2.402     12661     30.014    22.847
           3.603     18753     37.502    18.080
           4.804     28652     48.808    10.883
           6.005     42164     65.902    -0.000

      Averaged ODE under the yaw law: v_end = 3.074660086 km/s, di = 28.500000000 deg
      Direct optimisation (40 constant-yaw arcs, SLSQP): dV = 6.00538 km/s (converged), +0.001 % over Edelbaum

      Baselines:
        Sequential low-thrust (spiral, then plane change at GEO): 7.1120 km/s
        Impulsive Hohmann, optimal split (2.17 deg at perigee): 4.2670 km/s
        Low-thrust gravity-loss penalty vs. impulsive: 1.7383 km/s

      Cartesian propagation: m0 = 2000 kg, T = 2.0 N, Isp = 1800 s
        flight time 58.92 d (rocket-equation estimate 58.92 d), 532445 RHS evaluations
        propellant 576.8 kg (rocket equation 576.8 kg)
        final a = 42164.9 km (target 42164.2), e = 5.13e-03, i = 0.0087 deg

Most of the inclination is removed high up, where orbital speed is low: half
the delta-V is spent before 8 deg of plane change is done. Splitting the
work this way saves 1.1 km/s over spiralling first and turning at GEO. The
1.7 km/s premium over the impulsive transfer is the price of thrusting
continuously, paid back many times over by the 1800 s Isp: 577 kg of
propellant against roughly 1.5 t for a 320 s chemical stage on 4.27 km/s.
