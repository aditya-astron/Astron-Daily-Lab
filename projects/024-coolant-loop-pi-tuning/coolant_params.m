% Astron Built this AP
% Project 024 | PID and control systems | Coolant loop temperature control with transport delay: SIMC PI tuning
% Author: Aditya (aditya-astron)
% Date: 2026-11-09 IST
%
% Internal active thermal control loop (single-phase water/glycol), sized
% like an ISS moderate-temperature loop. A lumped coolant + cold-plate node T
% picks up the avionics load Q and rejects heat to an interface heat
% exchanger whose sink sits at Ts. A bypass valve sends a fraction u of the
% flow through the exchanger; with constant effectiveness
%
%     C dT/dt = Q - u G (T - Ts),        G = mdot cp eps
%
% The outlet RTD sits at the end of a pipe run, so it sees T delayed by the
% transport time theta = V_pipe / Vdot and filtered by its own lag tau_s.

function p = coolant_params()
  p.C      = 30e3;       % J/K, coolant inventory + cold plates
  p.mdot   = 0.06;       % kg/s
  p.cp     = 3560;       % J/(kg K), 50/50 propylene glycol
  p.rho    = 1040;       % kg/m^3
  p.eps    = 0.70;       % interface HX effectiveness
  p.Ts     = 4.0;        % degC, sink (low-temperature loop side)
  p.Tset   = 17.0;       % degC, outlet setpoint
  p.Q0     = 1500;       % W, nominal avionics load
  p.L_pipe = 6.0;        % m, valve to sensor
  p.d_pipe = 8e-3;       % m, inner diameter
  p.tau_s  = 4.0;        % s, RTD in thermowell

  p.G      = p.mdot * p.cp * p.eps;
  Vdot     = p.mdot / p.rho;
  p.theta  = pi * (p.d_pipe/2)^2 * p.L_pipe / Vdot;
  p.u0     = p.Q0 / (p.G * (p.Tset - p.Ts));   % trim valve position
end
