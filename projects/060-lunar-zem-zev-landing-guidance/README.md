# 060 - Lunar terminal descent with energy-optimal ZEM/ZEV guidance

## Problem

A lunar lander at high gate is 1.5 km up, 2 km short of the pad and 400 m
off to the side, moving at 45 m/s downrange and sinking at 25 m/s. It has a
45 kN throttleable engine (10 % floor, Isp 311 s) on a 12 t vehicle. Bring it
to rest on the pad, and choose the flight time so the engine never saturates
while spending as little propellant as the guidance law allows.

## Approach

Target-centred flat-Moon dynamics, valid over the last few kilometres:

    r'' = u + g,   g = (0, 0, -1.62) m/s^2,   m' = -|F| / (Isp g0),   u = F/m

**Guidance law.** Minimising int |u|^2 dt for a fixed time-to-go T with
r(T) = 0, v(T) = 0 gives a thrust acceleration that is linear in time.
Written with the zero-effort miss and velocity (engine off from now on):

    ZEM = -(r + v T + g T^2/2),   ZEV = -(v + g T)
    u   = 6 ZEM/T^2 - 2 ZEV/T
    J_E(T) = 12|ZEM|^2/T^3 - 12 ZEM.ZEV/T^2 + 4|ZEV|^2/T

Re-solving from any point on the arc returns the same line, so the law can
run every cycle in closed loop. Below 1 s to go the update stops (u grows
like 1/T^2) and the last profile a0 + a1 tau is flown open loop.

**Time-to-go.** J_E decreases monotonically with T, so a time penalty Gamma
is added: minimise Gamma T + J_E/2. Setting the derivative to zero gives

    (Gamma + |g|^2/2) T^4 - 2|v|^2 T^2 - 12 (r.v) T - 18|r|^2 = 0

with one sign change, hence exactly one positive root (Descartes).

**Thrust limit.** |u(t)| on the ideal arc is convex in t, so its peak is at
one end of the arc; the shortest feasible T is found by scanning and
bisecting on that peak against Tmax/m0 (conservative: the vehicle gets
lighter).

The simulator integrates position, velocity, mass and accumulated delta-v
with RK4 under a zero-order-hold guidance loop (10 Hz default) and clips
thrust to [0.1, 1.0] Tmax.

## Verification (`test_landing.py`)

- Vertical hand case: 1000 m above the pad, T = 40 s gives a0 = -2.13 m/s^2,
  a1 = 0.1875 m/s^3 and J_E = 292.476 m^2/s^3, all matched exactly.
- Closed-form kinematics of the profile hit r = 0, v = 0; J_E agrees with
  numerical quadrature; adding any perturbation orthogonal to {1, t} costs
  more energy.
- Re-solving at points along the arc reproduces the original line to 1e-10.
- Quartic root agrees with a bounded brute-force minimisation to 1e-5 s.
- With unlimited thrust the closed loop lands to round-off (1e-9 m, m/s)
  even at 1 Hz: each re-plan starts from the true state and the final
  segment is analytic, so sampling costs delta-v, not accuracy.
- Final mass equals m0 exp(-dv/ve) from the integrated delta-v to 1e-6 kg.
- At the feasible-T edge the limited engine never saturates and lands;
  at 0.7 of that T it saturates and arrives above 5 m/s.

## Results

| case | T [s] | dv [m/s] | propellant [kg] | peak throttle | landing |
|---|---|---|---|---|---|
| Gamma = 0.5 | 66.6 | 141.6 | 544.5 | 0.82 | exact |
| shortest feasible T (Gamma = 1.77) | 60.4 | 133.8 | 515.1 | 0.95 | exact |
| 0.8 x 66.6 s | 53.3 | 124.6 | 480.4 | 1.00 (86 cycles) | 10 m miss, 6.5 m/s, 8 m below pad |
| 3 x 66.6 s | 199.8 | 353.9 | 1314.7 | 0.53 | lands, but dips 4.4 m below pad level |

The delta-v minimum of the energy-optimal family (~127.6 m/s near T = 51 s)
needs 4.9 m/s^2 against the 3.75 m/s^2 available, so it is not flyable; past
the saturated region delta-v rises with T, which makes the feasibility edge
the best choice: 29 kg less propellant than an arbitrary Gamma = 0.5. Too
short a T saturates the engine late in the arc, where the fixed-time law
has no way to recover. Too long a T hovers against gravity and the
unconstrained arc can dip below the pad, which is the reason flight
versions add a glide-slope constraint.

## Run

    pip install numpy scipy
    python3 test_landing.py     # or: python3 -m pytest -q
    python3 study.py            # full study, see sample_output.txt
