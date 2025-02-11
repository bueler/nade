function [x, U, A] = generalbvp(p, q, f, a, b, alpha, beta, m)
% GENERALBVP  FD solution of the general linear 1D Dirichlet
% boundary value problem:
%   u''(x) + p(x) u'(x) + q(x) u(x) = f(x),  u(a) = alpha,  u(b) = beta
% Uses equally-spaced grid with m interior points, and centered
% finite differences.  Compare SECOND.
% Usage:
%   [x, U, A] = generalbvp(p, q, f, a, b, alpha, beta, m)

h = (b - a) / (m+1);           % grid spacing
x = (a:h:b)';                  % all grid points, including boundaries
xi = x(2:m+1);                 % interior points

A = sparse(m, m);              % use sparse matrix storage
FF = zeros(m, 1);
aa = -2 + h^2 * q(xi);         % diagonal entries
bb = 1 + h * p(xi) / 2;        % super-diagonal
cc = 1 - h * p(xi) / 2;        % sub-diagonal
for j = 1:m                    % assemble the jth row
    FF(j) = f(x(j+1));
    if j == 1
        A(1, [1, 2]) = [aa(1), bb(1)];
        FF(1) = FF(1) - alpha / h^2 + p(xi(1)) * alpha / (2 * h);
    elseif j == m
        A(m, [m-1, m]) = [cc(m), aa(m)];
        FF(m) = FF(m) - beta / h^2  - p(xi(m)) * beta  / (2 * h);
    else
        A(j, [j-1, j, j+1]) = [cc(j), aa(j), bb(j)];
    end
end
A = (1 / h^2) * A;

U = [alpha; A \ FF; beta];
