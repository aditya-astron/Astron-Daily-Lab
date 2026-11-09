% Astron Built this AP
% Project 024 | PID and control systems | Coolant loop temperature control with transport delay: SIMC PI tuning
% Author: Aditya (aditya-astron)
% Date: 2026-11-09 IST
%
% Linearize the coolant node about the trim point (Tset, u0, Q0):
%
%   d(dT)/dt = -(u0 G/C) dT - (G (Tset-Ts)/C) du + (1/C) dQ
%
% so tau = C/(u0 G), k_u = -(Tset-Ts)/u0 [K per unit valve], k_d = 1/(u0 G)
% [K/W]. The true plant is k e^{-theta s} / ((tau s + 1)(tau_s s + 1)).
% For tuning it is reduced to FOPDT with Skogestad's half rule: half of the
% smaller lag goes into the delay, half into the dominant lag.

function m = coolant_fopdt(p)
  m.tau    = p.C / (p.u0 * p.G);
  m.k      = -(p.Tset - p.Ts) / p.u0;
  m.kd     = 1 / (p.u0 * p.G);
  m.tau_s  = p.tau_s;
  m.theta0 = p.theta;                    % pure transport delay

  m.tau1   = m.tau + p.tau_s/2;          % half-rule FOPDT
  m.theta  = p.theta + p.tau_s/2;

  % frequency response of the full linear model (|k| because the
  % controller handles the sign as reverse action)
  m.mag   = @(w) abs(m.k) ./ (sqrt(1 + (w*m.tau).^2) .* sqrt(1 + (w*m.tau_s).^2));
  m.phase = @(w) -atan(w*m.tau) - atan(w*m.tau_s) - w*m.theta0;
end
