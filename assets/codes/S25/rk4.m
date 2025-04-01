function [tt,zz] = rk4(f,eta,t0,tf,N)
% RK4  Solve
%   u' = f(t,u),  u(t0) = eta
% for u(t) on the interval [t0,tf] with N steps of the classical
% O(k^4) Runge-Kutta method.
% Usage: [tt,zz] = rk4(f,eta,t0,tf,N)

dt = (tf - t0) / N;
tt = t0:dt:tf;        % row vector of times
eta = eta(:);         % force to be column vector
s = length(eta);
zz = zeros(s,N+1);    % jth column is U at t_{j-1}
zz(:,1) = eta;

for j = 1:N           % RK4 is (5.33) in LeVeque
    t = tt(j);
    Y1 = zz(:,j);
    f1 = f(t,Y1);
    Y2 = Y1 + (dt/2) * f1;
    f2 = f(t + dt/2,Y2);
    Y3 = Y1 + (dt/2) * f2;
    f3 = f(t + dt/2,Y3);
    Y4 = Y1 + dt * f3;
    f4 = f(t + dt,Y4);
    zz(:,j+1) = Y1 + (dt / 6) * (f1 + 2 * f2 + 2 * f3 + f4);
end
