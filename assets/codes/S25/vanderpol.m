% VANDERPOL  Use ODE45 to solve a famous nonlinear oscillation:
%   x'' - 5 (1-x^2) x' + x = 0
% This is written as a first order system with u(1) = x and u(2) = x'.
% Show the effects of adaptive time stepping.

f = @(t, u) [u(2); 5 * (1-u(1)^2) * u(2) - u(1)];

[tt, uu] = ode45(f, [0, 25], [2; 0]);

subplot(2,1,1)
plot(tt, uu(:,1), '.-')
ylabel('x(t)')

subplot(2,1,2)
plot(tt(2:end), diff(tt))
ylabel('\Delta t')
xlabel('t')
