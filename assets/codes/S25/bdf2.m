function [tt, zz] = bdf2(A, g, eta, t0, tf, N)
% BDF2  For a given square matrix A and a function g(t), solve
%   u'(t) = A u(t) + g(t),  u(t0) = eta
% for u(t) on the interval [t0,tf] with N equal timesteps of
% the multistep and implicit BDF2 method.  The first step uses
% the explicit trapezoid rule.  The combined scheme is
% O(k^2) (second-order) and A-stable.
% Usage: [tt, zz] = bdf2(A, g, eta, t0, tf, N)
% Example:
%   A = [0, 1; -4, 0];
%   g = @(t) zeros(2,1);
%   [tt, zz] = bdf2(A, g, [1; 0], 0, 5, 30);
%   plot(tt, zz(1,:), tt, cos(2*tt))
%   legend('BDF2', 'exact'), xlabel t

dt = (tf - t0) / N;
tt = t0:dt:tf;        % row vector of times
eta = eta(:);         % force to be column vector
s = length(eta);
zz = zeros(s, N+1);   % jth column is U at t_{j-1}
zz(:,1) = eta;

% first step: explicit trapezoid
f1 = A * zz(:,1) + g(tt(1));
zstar = zz(:,1) + dt * f1;
zz(:,2) = zz(:,1) + (dt/2) * (f1 + A * zstar + g(tt(2)));

% LU factor the relevant matrix in advance of solving for each step
% this replaces O(s^3) work at each step with O(s^2) work
M = 3 * eye(s,s) - 2 * dt * A;
[L, R, P] = lu(M);   % here P M = L R

% BDF2 multistep scheme, eqn (5.25) in LeVeque
for j = 2:N
    b = 4 * zz(:,j) - zz(:,j-1) + 2 * dt * g(tt(j+1));
    zz(:,j+1) = R \ (L \ P * b);   % solve  P M zz(:,j+1) = P b
end
