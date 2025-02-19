function [xx, yy, UU, A] = eastwind(m)
% EASTWIND  Solve the advection-diffusion equation
%   u_xx + u_yy - g(y) u_x = f(x,y)
% on the unit square.  There are zero Dirichlet boundary conditions
% on the north (y=1), south (y=0), and east (x=1) sides.  On the
% west side (x=0) there is a homogeneous Neumann condition u_x=0.
% The function g(y) gives advection which corresponds to a wind
% from the east.  The source term f(x,y) is a Gaussian which is
% concentrated around the point (0.8,0.4).
% Example with 2 visualizations:
%   [xx, yy, UU, A] = eastwind(80);
%   figure,  surf(xx, yy, UU),  view(2),  shading('interp'),  xlabel x,  ylabel y
%   figure,  mesh(xx, yy, UU),  xlabel x,  ylabel y
% For small m values, also:
%   figure,  spy(A)
% Compare:  POISSON

g = @(y) -40 * y * (1 - y);
f = @(x, y) - exp(-60 * ((x - 0.8)^2 + (y - 0.4)^2));

% grid suitable for *these* boundary conditions
h = 1.0 / (m + 1);
x = linspace(0, 1.0 - h, m + 1);
y = linspace(h, 1.0 - h, m);

% local-to-global grid index formula, correct for *this* grid
kk = @(i, j) (j - 1) * (m + 1) + i + 1;

% assemble system
N = m * (m + 1);
A = sparse(N, N);         % sparse matrix with no nonzero entries
F = zeros(N, 1);
for i = 0:m               % loop indices where there are unknowns U_ij
    for j = 1:m
        k = kk(i, j);     % assembing this *row*
        xi = x(i + 1);    % Matlab indexing starts at zero
        yj = y(j);
        F(k) = f(xi, yj);
        A(k, k) = -4.0;
        if i > 0          % west neighbor
            A(k, kk(i - 1, j)) = 1.0 + (h / 2) * g(yj);
        end
        if i < m          % east neighbor
            if i == 0
                A(k, kk(i + 1, j)) = 2.0; % from Neumann
            else
                A(k, kk(i + 1, j)) = 1.0 - (h / 2) * g(yj);
            end
        end
        if j > 1          % south neighbor
            A(k, kk(i, j - 1)) = 1.0;
        end
        if j < m          % north neighbor
            A(k, kk(i, j + 1)) = 1.0;
        end
    end
end
A = (1.0 / h^2) * A;

% solve
U = A \ F;

% put U on grid, including boundary values, for output
% note i and j loop ranges (local indices) are different from above;
% note Matlab indexing starts at zero!
UU = zeros(m + 2, m + 2);   % all grid points
for i = 0:m+1
    for j = 0:m+1
        if i == m+1 || j == 0 || j == m+1
            UU(i + 1, j + 1) = 0;   % zero Dirichlet values
        else
            UU(i + 1, j + 1) = U(kk(i, j));
        end
    end
end

% generate mesh for plotting
x = [x, 1.0];
y = [0.0, y, 1.0];
[xx, yy] = ndgrid(x, y);
