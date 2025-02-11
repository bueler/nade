function [x, U, A, FF] = neubvp(f, sig0, sig1, m)
% NEUBVP  Solve a Neumann-only boundary value problem.
% Note that f, sig0, sig1 need to satisfy a condition
% for a solution to exist.  Returns the linear system
% so that it can be edited.

h = 1.0 / (m+1);
x = (0:h:1)';

A = sparse(m+2, m+2);
FF = zeros(m+2, 1);
for j = 1:m+2
    if j == 1
        FF(1) = sig0 + (h/2) * f(x(1));
        A(1, [1, 2]) = [-1, 1] / h;
    elseif j == m+2
        FF(m+2) = - sig1 + (h/2) * f(x(m+2));
        A(m+2, [m+1, m+2]) = [1, -1] / h;
    else
        FF(j) = f(x(j));
        A(j, [j-1, j, j+1]) = [1, -2, 1] / h^2;
    end
end

U = A \ FF;
