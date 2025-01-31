function [A, x, U] = second(f, alpha, beta, m)
% SECONDA  Same as SECOND, but additionally returns A
% as the first output.

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
