function [tt,zz] = fe1(f,eta,t0,tf,N)
% FE1  Solve
%   u' = f(t,u),  u(t0) = eta
% for u(t) on the interval [t0,tf] with N steps of the forward
% Euler method.
% Usage: [tt,zz] = fe1(f,eta,t0,tf,N)

dt = (tf - t0) / N;
tt = t0:dt:tf;        % row vector of times
eta = eta(:);         % force to be column vector
s = length(eta);
zz = zeros(s,N+1);    % jth column is U at t_{j-1}
zz(:,1) = eta;

for j = 1:N           % forward Euler is (5.19) in LeVeque
    zz(:,j+1) = zz(:,j) + dt * f(tt(j),zz(:,j));
end
