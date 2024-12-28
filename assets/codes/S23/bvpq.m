function [x,U] = bvpq(m,xL,xR,q,f,alpha,beta)
% BVPQ  Solve the ODE BVP
%   u''(x) + q u(x) = f(x),  u(xL)=alpha,  u(xR)=beta
% using centered finite differences, for given input function f(x).
% Usage:
%   [x,u] = bvpq(m,xL,xR,q,f,alpha,beta)
% Example:  Solve an m=20 unknowns problem for which we know the
% exact solution.
%   >> f = @(x) (pi^2/4^2 + 1) * sin((pi/4)*x) - 1;
%   >> [x,U] = bvpq(20,0,2,-1,f,1,0);
%   >> uexact = @(x) 1 - sin((pi/4)*x);
%   >> plot(x,U,'o',x,uexact(x),'-')
%   >> xlabel x,  legend('numerical','exact')

% set up grid
h = (xR - xL) / (m+1);
x = xL:h:xR;           % length m+2

% assemble linear system  A U = F
A = (-2 + q * h^2) * eye(m);
for j = 1:m-1
    A(j,j+1) = 1;
    A(j+1,j) = 1;
end
A = (1/h^2) * A;
F = f(x(2:m+1))';    % evaluate f at interior grid points
F(1) = F(1) - alpha / h^2;
F(m) = F(m) - beta / h^2;

% solve the linear system
U = A \ F;           % numerical solution at interior points
U = [alpha U' beta]; % whole solution including boundary vals
