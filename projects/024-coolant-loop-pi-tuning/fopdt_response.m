% Astron Built this AP
% Project 024 | PID and control systems | Coolant loop temperature control with transport delay: SIMC PI tuning
% Author: Aditya (aditya-astron)
% Date: 2026-11-09 IST
%
% Frequency response of a plain FOPDT plant k e^{-theta s}/(tau s + 1),
% in the same mag/phase form as coolant_fopdt so the margin code is shared.

function m = fopdt_response(k, tau, theta)
  m.k = k; m.tau1 = tau; m.theta = theta;
  m.mag   = @(w) abs(k) ./ sqrt(1 + (w*tau).^2);
  m.phase = @(w) -atan(w*tau) - w*theta;
end
