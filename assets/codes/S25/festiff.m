% FESTIFF  Test forward Euler code FE1 on simple problem
%   u'(t) = - 25 u(t),  u(0) = 1,
% to compute u(0.8).  Tries three different step sizes
% and also plots the exact solution.
% Requires FE1.

f = @(t,u) -25 * u;

for N = [8 16 40]
    [tt, zz] = fe1(f, [1.0], 0, 0.8, N);
    plot(tt, zz, 'o-', 'linewidth', 0.5),  hold on
end
tfine = 0:0.001:0.8;
plot(tfine, exp(-25*tfine), 'k'),  hold off
axis([0 1 -1 1.2]),  xlabel t
legend('FE k=0.1', 'FE k=0.05', 'FE k=0.01', 'exact')