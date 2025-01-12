% LUMP  Simulate the temperature in a single lumped object.
% The model is Newton's law of cooling, an ordinary differential
% equation (ODE) initial value problem:
%   du/dt = k (Ta - u) + g(t),  u(0) = u0
% Here the heating source term g(t) is "on" over 1 < t < 2. (See
% values below.)  The numerical solution is extremely basic: the
% forward Euler method.  This code is NOT a good use of  Matlab.
% This simple ODE can be solved by hand, or very accurately by
% Matlab's built-in ODE45 solver.

u0 = 65.0;  % initial temperature
Ta = 65.0;  % ambient temperature
k = 0.2;    % conduction coefficient
g = @(t) (t > 1.0 && t < 2.0) * 80.0;  % source term
f = @(t, u) k * (Ta - u) + g(t);  % defines ODE

tf = 10.0;   % final time; total duration
NN = 200;   % number of steps

% Euler's method:
dt = tf / NN;  % time step
t = 0.0:dt:tf;
u = zeros(1, NN+1);
u(1) = u0;
for n = 1:NN
    u(n+1) = u(n) + dt * f(t(n), u(n));
end

% plot solution
plot(t, u, '.')
xlabel t,  ylabel('u(t) temperature'),  grid on
