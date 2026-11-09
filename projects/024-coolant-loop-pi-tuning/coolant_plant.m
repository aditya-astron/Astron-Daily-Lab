% Astron Built this AP
% Project 024 | PID and control systems | Coolant loop temperature control with transport delay: SIMC PI tuning
% Author: Aditya (aditya-astron)
% Date: 2026-11-09 IST
%
% Nonlinear coolant node plus RTD lag, packaged for simulate_loop.
% State x = [T; T_rtd_in]: T is the node temperature, the second state is
% the sensor lag driven by T. The pipe delay sits in the measurement line.
%
%   C dT/dt      = Q0 + dQ - u G (T - Ts)
%   tau_s dTm/dt = T - Tm
%
% Ordering note: delaying T then filtering equals filtering then delaying
% (both are LTI on T), so the lag can live in the state.

function pl = coolant_plant(p, dt)
  pl.f  = @(x, u, d) [ (p.Q0 + d - u*p.G*(x(1) - p.Ts)) / p.C ;
                       (x(1) - x(2)) / p.tau_s ];
  pl.h  = @(x) x(2);
  pl.x0 = [p.Tset; p.Tset];
  pl.nd = round(p.theta / dt);
end
