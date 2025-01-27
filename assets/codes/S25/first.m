% FIRST  First finite difference solution of steady heat equation in 1D.

alpha = 0.0;  beta = 0.0;      % boundary conditions
frhs = @(x) 1.0;               % solving  u''(x) = 1

m = 10;                        % resolution
h = 1.0 / (m+1);  x = 0:h:1;
A = zeros(m, m);  F = zeros(m, 1);
for j = 1:m                    % assemble the jth row
    F(j) = frhs(x(j+1));
    if j == 1
        A(1, [1, 2]) = [-2, 1] / h^2;
        F(1) = F(1) - alpha / h^2;
    elseif j == m
        A(m, [m-1, m]) = [1, -2] / h^2;
        F(m) = F(m) - beta / h^2;
    else
        A(j, [j-1, j, j+1]) = [1, -2, 1] / h^2;
    end
end

U = A \ F;
plot(x, [alpha, U', beta], 'o'),  xlabel x
