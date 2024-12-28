function [x,U] = bothneumann(m,f,sigma0,sigma1)
% BOTHNEUMANN  Attempt to solve the following ODE BVP, in which
% both ends of a rod have Neumann (derivative) boundary conditions:
%   u''(x) = f(x),  u'(0)=sigma0,  u'(1)=sigma1
% using centered finite differences, for given input function f(x).
% See section 2.13 of LeVeque for the method, and a discussion
% of why this is only an attempt.  THIS PROBLEM DOES NOT HAVE A
% UNIQUE SOLUTION.  Note the WARNING produced at run time:
%   warning: matrix singular to machine precision
% Example:  Attempt to solve m=20 test case problem.
%   >> f = @(x) 2 - 12 * x + 12 * x.^2;
%   >> [x,U] = bothneumann(20,f,0.0,0.0);
%   >> plot(x,U),  xlabel x
% This code is a small modification of STEADYHEAT, which solves
% a well-posed problem.

% input checking: is m = 1,2,3,...?
if m < 1 || floor(m) ~= m
    error('m must be an integer and at least 1')
end

% set up grid
h = 1.0 / (m+1);
x = 0:h:1;           % length m+2

% assemble linear system  A U = F
A = -2 * eye(m+2);
A(1,1:2) = [-h, h];
for j = 2:m+1
    A(j,[j-1,j+1]) = 1;
end
A(m+2,m+1:m+2) = [-h, h];
A = (1/h^2) * A;
% evaluate f at *all* grid points, then correct values
F = f(x)';
F(1) = sigma0 + (h/2) * f(x(1));
F(m+2) = sigma1 - (h/2) * f(x(m+2));

% solve the linear system
U = A \ F;           % numerical solution at interior points
