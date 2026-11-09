% Astron Built this AP
% Project 024 | PID and control systems | Coolant loop temperature control with transport delay: SIMC PI tuning
% Author: Aditya (aditya-astron)
% Date: 2026-11-09 IST
%
% SIMC PI rule (Skogestad 2003) for an FOPDT model, closed-loop time
% constant tau_c:
%
%   Kc = tau1 / (|k| (tau_c + theta)),   Ti = min(tau1, 4 (tau_c + theta))
%
% The min() keeps load-disturbance rejection reasonable on lag-dominant
% plants, where Ti = tau1 would leave a slow integrating tail.

function [Kc, Ti] = simc_pi(k, tau1, theta, tau_c)
  Kc = tau1 / (abs(k) * (tau_c + theta));
  Ti = min(tau1, 4 * (tau_c + theta));
end
