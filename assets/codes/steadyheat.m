function [x,U] = steadyheat(m,f,alpha,beta)
% STEADYHEAT  Solve the ODE BVP
%   u''(x) = f(x),  u(0)=alpha,  u(1)=beta
% using centered finite differences, for given input function f(x).
% (See section 2.4 of LeVeque for the method.)  Plottable outputs
% are the grid x and the numerical (approximate) solution U.  Both
% are length m+2 row vectors set-up for plotting.
% Usage:
%   [x, U] = steadyheat(m,f,alpha,beta)
% Example 1:  Solve an m=20 unknowns problem for which we know the
% exact solution.
%   >> f = @(x) 2 - 12 * x + 12 * x.^2;
%   >> [x,U] = steadyheat(20,f,0.0,0.0);
%   >> uexact = @(x) x.^2 .* (1 - x).^2;
%   >> plot(x,U,x,'o',uexact(x),'-')
%   >> xlabel x,  legend('numerical','exact')
% Example 2:  See VERIFYSH.

% input checking: is m = 1,2,3,...?
if m < 1 || floor(m) ~= m
    error('m must be an integer and at least 1')
end

% set up grid
h = 1.0 / (m+1);
x = 0:h:1;           % length m+2

% assemble linear system  A U = F
% WARNING: we should use sparse storage here
A = -2 * eye(m);
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
