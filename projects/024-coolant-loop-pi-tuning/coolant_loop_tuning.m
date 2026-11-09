% Astron Built this AP
% Project 024 | PID and control systems | Coolant loop temperature control with transport delay: SIMC PI tuning
% Author: Aditya (aditya-astron)
% Date: 2026-11-09 IST
%
% Tune the outlet-temperature PI loop of a coolant loop whose sensor sits
% behind a pipe delay, compare SIMC settings against Ziegler-Nichols, then
% fly the tunings on the nonlinear plant through an avionics load step, a
% setpoint change and a load that briefly exceeds the heat exchanger.

1;
function s = settle_time(t, y, band, t0)
  k = find(abs(y) > band, 1, 'last');
  if isempty(k), s = 0; else, s = t(k) - t0; end
end

p  = coolant_params();
m  = coolant_fopdt(p);
dt = 0.1;

printf('Coolant loop: C = %.0f kJ/K, G = %.1f W/K, Q0 = %.0f W, Tset = %.1f C, Ts = %.1f C\n', ...
       p.C/1e3, p.G, p.Q0, p.Tset, p.Ts);
printf('Trim valve u0 = %.3f, HX capacity at u = 1: %.0f W\n', p.u0, p.G*(p.Tset - p.Ts));
printf('Linear model: k = %.2f K/unit, kd = %.3f mK/W, tau = %.1f s, tau_s = %.1f s, pipe delay = %.2f s\n', ...
       m.k, 1e3*m.kd, m.tau, m.tau_s, m.theta0);
printf('Half-rule FOPDT: tau1 = %.1f s, theta = %.2f s  (tau1/theta = %.1f)\n', m.tau1, m.theta, m.tau1/m.theta);

[Ku, Pu, wu] = ultimate_gain(m);
printf('Ultimate gain Ku = %.3f /K, Pu = %.2f s\n\n', Ku, Pu);

names = {'SIMC tc=theta', 'SIMC tc=2theta', 'SIMC tc=4theta', 'Ziegler-Nichols'};
tc    = [1 2 4] * m.theta;
Kc = zeros(1, 4); Ti = zeros(1, 4);
for j = 1:3
  [Kc(j), Ti(j)] = simc_pi(m.k, m.tau1, m.theta, tc(j));
end
Kc(4) = 0.45*Ku; Ti(4) = Pu/1.2;

plant = coolant_plant(p, dt);
base  = struct('u0', p.u0, 'umin', 0, 'umax', 1, 'Tt', Inf, 'sgn', sign(m.k));

dQ = 300;
load_scn = struct('dt', dt, 'N', round(1500/dt), ...
                  'r', @(t) p.Tset, 'd', @(t) dQ*(t >= 30));
sp_scn   = struct('dt', dt, 'N', round(1500/dt), ...
                  'r', @(t) p.Tset - 1*(t >= 30), 'd', @(t) 0);

printf('Margins on the full model (exact delay + RTD lag); responses on the nonlinear plant\n');
printf('%-16s %7s %7s %6s %7s %6s %7s | %8s %8s %8s | %7s\n', 'tuning', 'Kc[/K]', 'Ti[s]', ...
       'GM', 'PM[deg]', 'Ms', 'DM[s]', 'peak[K]', 'IAE[Ks]', 'ts[s]', 'OS[%]');
for j = 1:4
  r = loop_margins(m, Kc(j), Ti(j));
  c = base; c.Kc = Kc(j); c.Ti = Ti(j);
  o = simulate_loop(plant, c, load_scn);
  err = o.x(:,1) - p.Tset;
  peak = max(err); iae = sum(abs(err))*dt;
  ts = settle_time(o.t, err, 0.02, 30);
  o2 = simulate_loop(plant, c, sp_scn);
  os = 100 * max(0, -(min(o2.x(:,1)) - (p.Tset - 1)));
  printf('%-16s %7.3f %7.1f %6.2f %7.1f %6.2f %7.2f | %8.3f %8.2f %8.0f | %7.1f\n', names{j}, ...
         Kc(j), Ti(j), r.GM, r.PM*180/pi, r.Ms, r.DM, peak, iae, ts, os);
end
printf('(load step +%d W at t = 30 s; ts = time to stay within 0.02 K; OS for a -1 K setpoint step)\n\n', dQ);

% Overload: the avionics load exceeds what the HX can reject at u = 1 for
% ten minutes. The valve pins open, T climbs, and an unprotected integrator
% keeps winding. When the load drops the stored integral holds the valve
% open and overcools the loop.
Qover = 2300;
ov = struct('dt', dt, 'N', round(2400/dt), 'r', @(t) p.Tset, ...
            'd', @(t) (Qover - p.Q0)*(t >= 60 & t < 660));
c = base; c.Kc = Kc(2); c.Ti = Ti(2);
printf('Overload %0.f W for 600 s (capacity %.0f W), SIMC tc=2theta\n', Qover, p.G*(p.Tset - p.Ts));
printf('%-22s %9s %11s %12s\n', 'integrator', 'Tmax[C]', 'Tmin[C]', 'recover[s]');
for Tt = [Inf, Ti(2)/2]
  c.Tt = Tt;
  o = simulate_loop(plant, c, ov);
  after = o.t >= 660;
  rec = settle_time(o.t(after), o.x(after,1) - p.Tset, 0.1, 660);
  if isinf(Tt), lbl = 'no anti-windup'; else, lbl = sprintf('back-calc Tt=%.0f s', Tt); end
  printf('%-22s %9.2f %11.2f %12.0f\n', lbl, max(o.x(:,1)), min(o.x(after,1)), rec);
end
