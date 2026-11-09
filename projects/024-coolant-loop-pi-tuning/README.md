# 024 - Coolant Loop Temperature Control with Transport Delay: SIMC PI Tuning

**Domain:** PID and control systems (active thermal control)

## Problem

An internal water/glycol coolant loop carries 1.5 kW of avionics heat to an
interface heat exchanger. A bypass valve sets the fraction `u` of the flow
that goes through the exchanger, and a PI loop holds the outlet at 17 C. The
RTD sits 6 m downstream, so it sees the coolant 5.2 s late and through its
own 4 s lag. Tune the loop, quantify how much robustness each tuning buys,
and check what happens when the load briefly exceeds the exchanger.

## Approach

Plant (lumped coolant + cold-plate node, constant HX effectiveness):

    C dT/dt = Q - u G (T - Ts),     G = mdot cp eps
    y(t) = RTD lag applied to T(t - theta),   theta = V_pipe / Vdot

Linearized at trim `u0 = Q0 / (G (Tset - Ts))`:

    tau = C/(u0 G) = 260 s,   k = -(Tset - Ts)/u0 = -16.8 K/unit,   kd = 1/(u0 G)

The 4 s sensor lag is folded into an FOPDT model with Skogestad's half rule
(`tau1 = tau + tau_s/2`, `theta = theta_pipe + tau_s/2`), then tuned with SIMC:

    Kc = tau1 / (|k| (tau_c + theta)),   Ti = min(tau1, 4 (tau_c + theta))

Margins are evaluated on the full model with the exact delay (phases summed
analytically, no Pade): GM, PM, delay margin `PM/wc` and peak sensitivity
`Ms`. Ziegler-Nichols PI from the ultimate gain is the comparison. The loop
is then simulated on the nonlinear plant (RK4, ZOH valve, sampled delay
line, 0.1 s) with back-calculation anti-windup available.

## Verification (`test_coolant_loop.m`, 24 checks)

- SIMC `tau_c = theta` on `tau1 <= 8 theta` gives `L = e^{-theta s}/(2 theta s)`:
  GM = pi, PM = 61.35 deg, DM = (pi - 1) theta, Ms = 1.59 (Skogestad 2003).
- Ultimate gain: P control at 0.97/1.00/1.03 Ku decays / holds / grows in
  simulation, oscillating at the computed Pu to 0.04%.
- Delay margin: the PI loop is stable at 0.95 DM extra delay and unstable at 1.05 DM.
- Linearization (tau, k, kd) matches central differences; pipe delay = V/Vdot.
- Delay line shifts the sensor state by exactly `nd` samples.
- Integral-error identity on the nonlinear plant: with no saturation,
  `int e dt = Ti/Kc (u1 - u0)` where `u1` comes from the new heat balance
  (matches to 1e-6); zero steady-state error.
- Back-calculation cuts post-overload overcooling by more than 10x.

## How to run

```
octave --no-gui -q coolant_loop_tuning.m   # tuning study
octave --no-gui -q test_coolant_loop.m     # checks (~1.5 min)
```

Tested with GNU Octave 8 (no packages needed).

## Sample output

```
Coolant loop: C = 30 kJ/K, G = 149.5 W/K, Q0 = 1500 W, Tset = 17.0 C, Ts = 4.0 C
Trim valve u0 = 0.772, HX capacity at u = 1: 1944 W
Linear model: k = -16.85 K/unit, kd = 8.667 mK/W, tau = 260.0 s, tau_s = 4.0 s, pipe delay = 5.23 s
Half-rule FOPDT: tau1 = 262.0 s, theta = 7.23 s  (tau1/theta = 36.2)
Ultimate gain Ku = 3.512 /K, Pu = 34.26 s

Margins on the full model (exact delay + RTD lag); responses on the nonlinear plant
tuning            Kc[/K]   Ti[s]     GM PM[deg]     Ms   DM[s] |  peak[K]  IAE[Ks]    ts[s] |   OS[%]
SIMC tc=theta      1.076    57.8   2.95    43.0   1.85   10.84 |    0.148     8.29      107 |    52.2
SIMC tc=2theta     0.717    86.7   4.59    56.2   1.46   20.93 |    0.185    18.66      199 |    35.6
SIMC tc=4theta     0.430   144.6   7.86    69.1   1.24   42.66 |    0.264    51.84      421 |    13.1
Ziegler-Nichols    1.580    28.5   1.76    21.0   3.39    3.64 |    0.123     4.09       64 |    69.5
(load step +300 W at t = 30 s; ts = time to stay within 0.02 K; OS for a -1 K setpoint step)

Overload 2300 W for 600 s (capacity 1944 W), SIMC tc=2theta
integrator               Tmax[C]     Tmin[C]   recover[s]
no anti-windup             19.27       14.30          842
back-calc Tt=43 s          19.27       17.00          176
```

## Notes

- The plant is lag-dominant (`tau1/theta = 36`), so SIMC cuts `Ti` to
  `4 (tau_c + theta)`. That is what keeps the load-step IAE small, and it is
  also why setpoint steps overshoot; a setpoint weight or prefilter fixes
  that without touching disturbance rejection.
- `tau_c = 2 theta` is the practical pick: Ms ~ 1.5 and 21 s of delay margin,
  four times the 5.2 s pipe delay, so a pump running at half flow (which
  doubles that delay) stays well inside
  it, while the loop still rejects a 300 W step with under 0.2 K excursion.
- Ziegler-Nichols reacts fastest but at 21 deg PM and Ms = 3.4 it is one
  pump degradation away from oscillating.
- With no anti-windup the integrator keeps winding while the valve is pinned
  open; when the load drops it overcools the loop by 2.7 K and takes 14 min
  to recover.
- The IAE for the load step equals the integral-error identity value
  `Ti/Kc (u1 - u0)` whenever the response does not cross the setpoint,
  which is why it scales with `Ti/Kc`, not with the tuning's speed alone.
