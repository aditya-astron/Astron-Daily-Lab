% Astron Built this AP
% Project 024 | PID and control systems | Coolant loop temperature control with transport delay: SIMC PI tuning
% Author: Aditya (aditya-astron)
% Date: 2026-11-09 IST
%
% Ultimate (critical) proportional gain: the frequency where the plant
% phase reaches -pi and the gain that puts |L| = 1 there. Ziegler-Nichols
% PI from these: Kc = 0.45 Ku, Ti = Pu / 1.2.

function [Ku, Pu, wu] = ultimate_gain(m)
  wu = fzero(@(w) m.phase(w) + pi, [1e-9, pi / m.theta + 10]);
  Ku = 1 / m.mag(wu);
  Pu = 2*pi / wu;
end
