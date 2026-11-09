% Astron Built this AP
% Project 024 | PID and control systems | Coolant loop temperature control with transport delay: SIMC PI tuning
% Author: Aditya (aditya-astron)
% Date: 2026-11-09 IST
%
% Sampled PI loop around a continuous plant with a measurement delay.
%
% plant: f(x,u,d) state derivative, h(x) measured quantity before the delay,
%        x0 initial state, nd delay in samples (theta = nd*dt).
% ctrl:  Kc, Ti, u0 (bias), umin/umax, Tt (back-calculation tracking time,
%        Inf = no anti-windup), sgn = sign of the plant gain.
% scen:  dt, N, r(t) setpoint, d(t) disturbance.
%
% Error is formed as e = sgn (r - y) so a negative-gain plant (opening the
% HX valve cools) gets reverse action with a positive Kc. The integrator is
% held in valve units:
%
%   u  = u0 + Kc e + I,   I <- I + dt (Kc/Ti e + (u_sat - u)/Tt)
%
% The plant is advanced with RK4 under a zero-order-hold input.

function out = simulate_loop(plant, ctrl, scen)
  dt = scen.dt; N = scen.N;
  x = plant.x0;
  ybuf = repmat(plant.h(x), plant.nd + 1, 1);   % circular delay line
  head = 1;
  I = 0;

  out.t = (0:N-1)' * dt;
  [out.y, out.u, out.e, out.r, out.d] = deal(zeros(N, 1));
  out.x = zeros(N, numel(x));

  for n = 1:N
    t = out.t(n);
    y = ybuf(head);                 % h(x(t - theta))
    r = scen.r(t); d = scen.d(t);
    e = ctrl.sgn * (r - y);

    u_raw = ctrl.u0 + ctrl.Kc*e + I;
    u = min(max(u_raw, ctrl.umin), ctrl.umax);
    I = I + dt * (ctrl.Kc/ctrl.Ti * e + (u - u_raw)/ctrl.Tt);

    out.y(n) = y; out.u(n) = u; out.e(n) = e; out.r(n) = r; out.d(n) = d;
    out.x(n, :) = x(:)';

    k1 = plant.f(x, u, d);
    k2 = plant.f(x + dt/2*k1, u, d);
    k3 = plant.f(x + dt/2*k2, u, d);
    k4 = plant.f(x + dt*k3, u, d);
    x = x + dt/6*(k1 + 2*k2 + 2*k3 + k4);

    ybuf(head) = plant.h(x);        % sample taken at t + dt
    head = mod(head, plant.nd + 1) + 1;
  end
  out.I_end = I;
end
