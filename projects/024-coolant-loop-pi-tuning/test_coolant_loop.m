% Astron Built this AP
% Project 024 | PID and control systems | Coolant loop temperature control with transport delay: SIMC PI tuning
% Author: Aditya (aditya-astron)
% Date: 2026-11-09 IST
%
% Verification: analytic SIMC margins, margins vs. time-domain stability
% boundaries, linearization, delay line, and the PI integral-error identity.

1;
global npass nfail
npass = 0; nfail = 0;

function check(name, ok, detail)
  global npass nfail
  if ok, npass++; s = 'PASS'; else, nfail++; s = 'FAIL'; end
  printf('[%s] %-58s %s\n', s, name, detail);
end

function g = growth(y, N)
  % peak |y| over the last fifth of the record vs. the fifth before it
  a = max(abs(y(round(0.6*N)+1:round(0.8*N))));
  b = max(abs(y(round(0.8*N)+1:N)));
  g = b / a;
end

% 1. SIMC with tau_c = theta on tau1 <= 8 theta cancels the lag (Ti = tau1)
%    and leaves L = e^{-theta s}/(2 theta s): GM = pi, PM = pi/2 - 1/2 rad,
%    DM = (pi - 1) theta, Ms = 1.59 (Skogestad 2003, Table 1).
th = 2.0; tau = 5*th; k = 3.0;
m = fopdt_response(k, tau, th);
[Kc, Ti] = simc_pi(k, tau, th, th);
r = loop_margins(m, Kc, Ti);
check('SIMC tc=theta: Ti = tau1', abs(Ti - tau) < 1e-12, sprintf('Ti = %.3f', Ti));
check('GM = pi', abs(r.GM - pi) < 1e-9, sprintf('%.12f', r.GM));
check('PM = pi/2 - 0.5 rad (61.35 deg)', abs(r.PM - (pi/2 - 0.5)) < 1e-9, sprintf('%.6f deg', r.PM*180/pi));
check('DM = (pi - 1) theta', abs(r.DM - (pi - 1)*th) < 1e-8, sprintf('%.6f s', r.DM));
check('Ms = 1.59', abs(r.Ms - 1.59) < 0.005, sprintf('%.4f', r.Ms));

% 2. Ultimate gain: phase condition holds, and P control at Ku sits on the
%    stability boundary in simulation, oscillating at period Pu.
[Ku, Pu, wu] = ultimate_gain(m);
check('atan(wu tau) + wu theta = pi', abs(atan(wu*tau) + wu*th - pi) < 1e-12, sprintf('wu = %.6f rad/s', wu));

% Sampling with ZOH adds about dt/2 of effective delay, so the "neutral"
% case drifts slightly unstable as dt grows; 2 ms keeps that under 2%.
dt = 2e-3; N = round(400/dt);
pl = fopdt_plant(k, tau, th, dt);
scn = struct('dt', dt, 'N', N, 'r', @(t) 1*(t < 0.5), 'd', @(t) 0);
c = struct('u0', 0, 'umin', -Inf, 'umax', Inf, 'Tt', Inf, 'sgn', 1, 'Ti', Inf);
g = zeros(1, 3); f = [0.97 1.00 1.03];
for i = 1:3
  c.Kc = f(i)*Ku;
  o = simulate_loop(pl, c, scn);
  g(i) = growth(o.y, N);
  if i == 2
    y = o.y(round(0.5*N):end);
    z = find(y(1:end-1) < 0 & y(2:end) >= 0);
    Psim = mean(diff(z))*dt;
  end
end
check('P control at 0.97 Ku decays', g(1) < 0.99, sprintf('growth %.4f', g(1)));
check('P control at 1.00 Ku neutral', abs(g(2) - 1) < 0.02, sprintf('growth %.4f', g(2)));
check('P control at 1.03 Ku grows', g(3) > 1.01, sprintf('growth %.4f', g(3)));
check('oscillation period = Pu', abs(Psim - Pu)/Pu < 2e-3, sprintf('%.4f vs %.4f s', Psim, Pu));

% 3. Delay margin from loop_margins vs. simulation of the PI loop with the
%    delay stretched just under and just over theta + DM.
r = loop_margins(m, Kc, Ti);
c = struct('u0', 0, 'umin', -Inf, 'umax', Inf, 'Tt', Inf, 'sgn', 1, 'Kc', Kc, 'Ti', Ti);
N = round(600/dt);
scn.N = N;
for fdm = [0.95 1.05]
  o = simulate_loop(fopdt_plant(k, tau, th + fdm*r.DM, dt), c, scn);
  gg = growth(o.y, N);
  if fdm < 1
    check(sprintf('PI stable with %.2f x delay margin', fdm), gg < 1, sprintf('growth %.4f', gg));
  else
    check(sprintf('PI unstable with %.2f x delay margin', fdm), gg > 1, sprintf('growth %.4f', gg));
  end
end

% 4. Coolant linearization against central differences of the nonlinear
%    node equation at trim.
p  = coolant_params();
mc = coolant_fopdt(p);
fT = @(T, u, Q) (Q - u*p.G*(T - p.Ts)) / p.C;
h = 1e-4;
dfdT = (fT(p.Tset+h, p.u0, p.Q0) - fT(p.Tset-h, p.u0, p.Q0)) / (2*h);
dfdu = (fT(p.Tset, p.u0+h, p.Q0) - fT(p.Tset, p.u0-h, p.Q0)) / (2*h);
dfdQ = (fT(p.Tset, p.u0, p.Q0+1) - fT(p.Tset, p.u0, p.Q0-1)) / 2;
check('trim is an equilibrium', abs(fT(p.Tset, p.u0, p.Q0)) < 1e-15, '');
check('tau = -1/(df/dT)', abs(-1/dfdT - mc.tau)/mc.tau < 1e-9, sprintf('%.4f s', mc.tau));
check('k = (df/du)/(-df/dT)', abs(dfdu/(-dfdT) - mc.k)/abs(mc.k) < 1e-9, sprintf('%.4f K', mc.k));
check('kd = (df/dQ)/(-df/dT)', abs(dfdQ/(-dfdT) - mc.kd)/mc.kd < 1e-9, sprintf('%.6f K/W', mc.kd));
check('pipe delay = V/Vdot', abs(mc.theta0 - pi*0.004^2*6/(0.06/1040)) < 1e-12, sprintf('%.4f s', mc.theta0));

% 5. Delay line: open loop with a constant valve offset, measured output
%    equals the sensor state shifted by exactly nd samples.
dt = 0.1;
pl = coolant_plant(p, dt);
c = struct('u0', p.u0 - 0.1, 'umin', 0, 'umax', 1, 'Tt', Inf, 'sgn', -1, 'Kc', 0, 'Ti', Inf);
o = simulate_loop(pl, c, struct('dt', dt, 'N', 3000, 'r', @(t) p.Tset, 'd', @(t) 0));
nd = pl.nd;
check('measurement = sensor state delayed nd samples', ...
      max(abs(o.y(nd+1:end) - o.x(1:end-nd, 2))) == 0, sprintf('nd = %d', nd));

% 6. Integral-error identity on the nonlinear plant. Without saturation the
%    final valve move is carried by the integrator alone, Du = (Kc/Ti) int e,
%    and Du is fixed by the new heat balance: u1 = (Q0 + dQ)/(G (Tset - Ts)).
%    So int (T_meas - Tset) dt = Ti/Kc * (u1 - u0) regardless of dynamics.
[Kc, Ti] = simc_pi(mc.k, mc.tau1, mc.theta, 2*mc.theta);
dQ = 300;
c = struct('u0', p.u0, 'umin', 0, 'umax', 1, 'Tt', Inf, 'sgn', -1, 'Kc', Kc, 'Ti', Ti);
o = simulate_loop(pl, c, struct('dt', dt, 'N', round(4000/dt), 'r', @(t) p.Tset, 'd', @(t) dQ*(t >= 30)));
u1 = (p.Q0 + dQ) / (p.G*(p.Tset - p.Ts));
ie_sim = sum(o.e)*dt;
ie_th  = Ti/Kc * (u1 - p.u0);
check('int e dt = Ti/Kc (u1 - u0)', abs(ie_sim - ie_th)/ie_th < 1e-6, sprintf('%.6f vs %.6f K s', ie_sim, ie_th));
check('no steady-state error after load step', abs(o.x(end,1) - p.Tset) < 1e-6, ...
      sprintf('%.2e K', o.x(end,1) - p.Tset));
check('valve settles at heat-balance position', abs(o.u(end) - u1) < 1e-6, sprintf('%.6f', o.u(end)));
check('valve never saturated', all(o.u > 0 & o.u < 1), '');

% 7. Back-calculation anti-windup bounds the overcool after an overload.
ov = struct('dt', dt, 'N', round(2400/dt), 'r', @(t) p.Tset, 'd', @(t) 800*(t >= 60 & t < 660));
c.Tt = Inf;    o1 = simulate_loop(pl, c, ov);
c.Tt = Ti/2;   o2 = simulate_loop(pl, c, ov);
a = o1.t >= 660;
under1 = p.Tset - min(o1.x(a,1)); under2 = max(0, p.Tset - min(o2.x(a,1)));
check('valve saturates during overload', max(o1.u) == 1, '');
check('anti-windup cuts overcool by > 10x', under2 < under1/10, sprintf('%.3f vs %.3f K', under2, under1));

printf('\n%d passed, %d failed\n', npass, nfail);
if nfail > 0, exit(1); end
