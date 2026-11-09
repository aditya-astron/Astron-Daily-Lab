% Astron Built this AP
% Project 024 | PID and control systems | Coolant loop temperature control with transport delay: SIMC PI tuning
% Author: Aditya (aditya-astron)
% Date: 2026-11-09 IST
%
% Stability margins of L = C_PI(jw) G(jw) with the exact (irrational) delay.
%
%   C_PI = Kc (1 + 1/(jw Ti)):  |C| = Kc sqrt(1 + 1/(w Ti)^2),
%                                arg C = atan(w Ti) - pi/2
%
% Phases are summed analytically, so no unwrapping is needed.
%   gain crossover wc: |L| = 1,       PM = pi + arg L(wc)
%   phase crossover w180: arg L = -pi, GM = 1/|L(w180)|
%   delay margin DM = PM / wc (extra dead time the loop tolerates)
%   Ms = max |1/(1+L)|, peak sensitivity

function r = loop_margins(m, Kc, Ti)
  magL = @(w) Kc * sqrt(1 + 1./(w*Ti).^2) .* m.mag(w);
  phL  = @(w) atan(w*Ti) - pi/2 + m.phase(w);

  w = logspace(-6, 2, 20000);

  i = find(magL(w) < 1, 1);
  r.wc = fzero(@(x) log(magL(x)), [w(i-1), w(i)]);
  r.PM = pi + phL(r.wc);

  i = find(phL(w) < -pi, 1);
  r.w180 = fzero(@(x) phL(x) + pi, [w(i-1), w(i)]);
  r.GM = 1 / magL(r.w180);

  r.DM = r.PM / r.wc;

  S = @(x) 1 ./ abs(1 + magL(x) .* exp(1i*phL(x)));
  [~, i] = max(S(w));
  lo = w(max(i-1, 1)); hi = w(min(i+1, numel(w)));
  wpk = fminbnd(@(x) -S(x), lo, hi, optimset('TolX', 1e-12));
  r.Ms = S(wpk);
end
