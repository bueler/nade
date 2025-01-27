function [x, U] = second(f, alpha, beta, m)
% SECOND  A better finite difference solution of a steady heat equation
% problem in 1D:
%   u''(x) = f(x),  u(0) = alpha,  u(1) = beta
% Uses equally-spaced grid with m interior points, and centered
% finite differences.  Compare FIRST.
% Usage:
%   [x, U] = second(f, alpha, beta, m)
% Example (actually *too easy*; compare THIRD):
%   >> f = @(x) 1;
%   >> uexact = @(x) (x/2) .* (x-1);
%   >> m = 10;
%   >> [x, U] = second(f, 0.0, 0.0, m);
%   >> plot(x, U, 'o', x, uexact(x))
%   >> xlabel x,  legend('U_j', 'u(x)')
%   >> max(abs(U - uexact(x)))            % numerical error

h = 1.0 / (m+1);
x = (0:h:1)';

A = sparse(m, m);              % use sparse matrix storage
FF = zeros(m, 1);
for j = 1:m                    % assemble the jth row
    FF(j) = f(x(j+1));
    if j == 1
        A(1, [1, 2]) = [-2, 1] / h^2;
        FF(1) = FF(1) - alpha / h^2;
    elseif j == m
        A(m, [m-1, m]) = [1, -2] / h^2;
        FF(m) = FF(m) - beta / h^2;
    else
        A(j, [j-1, j, j+1]) = [1, -2, 1] / h^2;
    end
end

U = [alpha; A \ FF; beta];
