function [x,y,A,b] = passemble(m,g)
% ASSEMBLE  Assemble matrix A and right hand side b, in
% linear system  A u = b,  for the Poisson equation
%   u_xx + u_yy = g(x,y)
% on the unit square, with zero Dirichlet boundary conditions.
% Usage:  [x,y,A,b] = passemble(m,f)

% grid and indexing
h = 1.0 / (m+1);
x = linspace(h,1.0-h,m);  y = x;
kk = @(i,j) (j-1) * m + i;     % local-to-global indexing

% assemble linearization; note  A u ~ u_xx + u_yy
N = m^2;  A = sparse(N,N);
b = zeros(N,1);
for i = 1:m
    for j = 1:m
        k = kk(i,j);
        A(k,k) = -4.0;
        if i > 1,    A(k,kk(i-1,j)) = 1.0;    end
        if i < m,    A(k,kk(i+1,j)) = 1.0;    end
        if j > 1,    A(k,kk(i,j-1)) = 1.0;    end
        if j < m,    A(k,kk(i,j+1)) = 1.0;    end
        b(k) = g(x(i),y(j));
    end
end
A = (1.0 / h^2) * A;
