function [xx,yy,UU,A] = poisson(m,f)
% POISSON  Solve the Poisson equation
%   u_xx + u_yy = f(x,y)
% on the unit square, with zero Dirichlet boundary conditions.
% Usage:   [x,y,UU,A] = poisson(m,f)
% Example with general right-hand-side and surface plot:
%   >> f = @(x,y) exp(-30.0*((x-0.5).^2 + (y-0.75).^2));
%   >> [x,y,U] = poisson(100,f);
%   >> surf(x,y,U), xlabel x, ylabel y
% Example of verification case:  POISSONCONV
% Example showing sparsity pattern of A:
%   >> [x,y,U,A] = poisson(5);  spy(A)

% grid
h = 1.0 / (m+1);
x = linspace(h,1.0-h,m);  y = x;
[xx,yy] = ndgrid(x,y);

% assemble system
kk = @(i,j) (j-1) * m + i;    % local-to-global grid index formula
N = m^2;  A = sparse(N,N);  F = zeros(N,1);
for i = 1:m
    for j = 1:m
       k = kk(i,j);
       F(k) = f(x(i),y(j));
       A(k,k) = -4.0;
       if i > 1,    A(k,kk(i-1,j)) = 1.0;    end
       if i < m,    A(k,kk(i+1,j)) = 1.0;    end
       if j > 1,    A(k,kk(i,j-1)) = 1.0;    end
       if j < m,    A(k,kk(i,j+1)) = 1.0;    end
    end
end
A = (1.0 / h^2) * A;
U = A \ F;           % solve

% put U on grid for output
UU = zeros(m,m);
for i = 1:m
    for j = 1:m
        UU(i,j) = U(kk(i,j));
    end
end
