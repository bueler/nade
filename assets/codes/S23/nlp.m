function [x,y,U,err] = nlp(m,gamma)
% NLP  Solve the nonlinear Poisson equation
%   u_xx + u_yy + gamma u^3 = g(x,y)
% on the unit square, with zero Dirichlet boundary conditions.
% Uses manufactured solution to compute error.
% Calls PASSEMBLE.
% Usage:   [x,y,UU,err] = nlp(m,gamma)
% Example:
%   >> [x,y,UU,err] = nlp(10,10.0);
%   >> surf(x,y,UU)

uexact = @(x,y) sin(pi * x) .* sin(2 * pi * y);
g = @(x,y) - 5 * pi^2 * uexact(x,y) + gamma * uexact(x,y).^3;

[x,y,A,b] = passemble(m,g);
h = 1.0 / (m+1);
N = m^2;

% Newton iteration
U = zeros(N,1);
for knewt = 1:15
    FF = A * U + gamma * U.^3 - b;
    rnorm = norm(FF);
    printf('    %2d:  residual norm %.4e\n',knewt-1,rnorm)
    if knewt == 1,  rnorm0 = rnorm;  end
    if (rnorm / rnorm0) < 1.0e-9,  break;  end
    JJ = A + 3 * gamma * spdiags(U.^2,0,N,N);
    s = - JJ \ FF;                % Newton step
    U = U + s;
end

[xx,yy] = ndgrid(x,y);
Uex = reshape(uexact(xx,yy),N,1);
err = norm(U-Uex,Inf);
printf('on m=%d grid, h = %.5f:  |U-Uexact|_inf = %.3e\n',...
       m,h,err);
U = reshape(U,m,m);
