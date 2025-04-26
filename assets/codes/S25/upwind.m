function [x, u] = upwind(a, tf, h)
% UPWIND  Solve the constant-coefficient, a >= 0 advection equation
%   u_t + a u_x = 0
% on [0,1] for u(t,x).  Boundary condition is u(t,0) = 0.  The
% initial condition u(0,x) is equal to the characteristic function
% of (0.1,0.3).  This function returns the grid x and the final-time
% solution u(tf,x).  (This function does no plotting.)
% Usage:  [x, u] = upwind(a, tf, h)
% Example:
%   >> [x, u]=upwind(1, 0.5, 0.01);
%   >> plot(x, eta, x, u)

m = ceil(1.0 / h);
h = 1.0 / m;

if a < 0, error('only written for a >= 0'), end
if a == 0
    k = h;            % irrelevant ... no motion
else
    k = 0.8 * h / a;  % 80% of CFL time step
end
NN = ceil(tf / k);
k = tf / NN;          % guarantees that  k NN = tf

x = 0.0:h:1.0;
eta = zeros(size(x));
mark = (x > 0.1) & (x < 0.3);
eta(mark) = 1.0;

fprintf('doing N=%d time steps to tf=%.4f ...\n', NN, tf)
r = a * k / h;
u = eta;
u(1) = 0.0;   % unaltered boundary condition
for n = 1:NN
    u(2:end) = u(2:end) - r * (u(2:end) - u(1:end-1));
end
