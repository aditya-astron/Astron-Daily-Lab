% Astron Built this AP
% Project 024 | PID and control systems | Coolant loop temperature control with transport delay: SIMC PI tuning
% Author: Aditya (aditya-astron)
% Date: 2026-11-09 IST
%
% Linear FOPDT plant in deviation variables for simulate_loop:
% tau dx/dt = -x + k u + kd d, measured x delayed by theta.

function pl = fopdt_plant(k, tau, theta, dt, kd)
  if nargin < 5, kd = 1; end
  pl.f  = @(x, u, d) (-x + k*u + kd*d) / tau;
  pl.h  = @(x) x;
  pl.x0 = 0;
  pl.nd = round(theta / dt);
end
