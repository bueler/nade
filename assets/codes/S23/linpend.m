function [x,U] = linpend(m,a,b,alpha,beta)
% LINPEND  Solve the ODE BVP
%   u''(x) + u(x) = 0,  u(a)=alpha,  u(b)=beta
% using centered finite differences.  Example: TESTLINPEND

% set up grid
h = (b - a) / (m+1);
x = a:h:b;            % length m+2

% assemble linear system  A U = F  using sparse storage
A = -2 * speye(m);
for j = 1:m-1
    A(j,j+1) = 1;
    A(j+1,j) = 1;
end
A = (1/h^2) * A + speye(m);
F = zeros(m,1);
F(1) = - alpha / h^2;
F(m) = - beta / h^2;

% solve the linear system
U = A \ F;           % numerical solution at interior points
U = [alpha U' beta]; % whole solution including boundary vals
