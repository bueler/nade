function [xx,yy,UU] = heat2d(m,fsource)
% HEAT2D  Solve Poisson equation on the unit square [0,1]^2:
%   u_xx + u_yy = f(x,y)
% with Dirichlet conditions u=0 on the x=0, x=1, and y=1 sides, but
% u_y=0 on the y=0 side.  The grid has local (left fig.) indices,
% while the unknowns U_k = * have global (right fig.) indices as
% follows in the m = 3 case.  Note the Dirichlet boundary is
% indicated by a zero:
%
%   y=1  j=4  0   0   0   0   0           0   0   0   0   0
%        j=3  0   *   *   *   0           0   10  11  12  0
%        j=2  0   *   *   *   0           0   7   8   9   0
%        j=1  0   *   *   *   0           0   4   5   6   0
%   y=0  j=0  0   *   *   *   0           0   1   2   3   0
%            i=0 i=1 i=2 i=3 i=4
%            x=0             x=1              k values
%
% Note h = 1/(m+1) and (x_i,y_j) = (i h,j h) is grid location.
% Usage:
%    >> [xx,yy,UU] = heat2d(m,fsource)
% Example:  See TESTHEAT2D.

h = 1 / (m+1);
n = m * (m+1);                       % number of unknowns
A = sparse(n,n);
F = zeros(n,1);
%KK = zeros(size(xx));               % for debugging: store global indices
for j = 0:m
    for i = 1:m
        k = LG(m,i,j);               % at (x_i,y_j) --> U_k is unknown here
        %KK(i+1,j+1) = k;            % for debugging
        F(k) = fsource(i*h,j*h);
        % fill kth row of A
        A(k,k) = -4;
        if i > 1
            A(k,LG(m,i-1,j)) = 1;    % west neighbor
        end
        if i < m
            A(k,LG(m,i+1,j)) = 1;    % east neighbor
        end
        if j > 0
            A(k,LG(m,i,j-1)) = 1;    % south neighbor
        end
        if j == 0
            A(k,LG(m,i,j+1)) = 2;    % north neighbor (Neumann case)
        elseif j < m
            A(k,LG(m,i,j+1)) = 1;    % north neighbor (general)
        end
    end
end
A = (1/h^2) * A;
U = A \ F;                           % solution as a column vector

% for debugging: to check matrix structure
%spy(A)
%full(A)
% for debugging: to check grid ordering, uncomment *all* "KK" lines
%KK, xx, yy

% insert solution into grid positions, in a form suitable
% for plotting
[xx,yy] = ndgrid(0:h:1,0:h:1);       % grid locations
UU = zeros(size(xx));                % includes Dirichlet values
UU(2:m+1,1:m+1) = reshape(U,m,m+1);  % insert
end % heat2d

    function k = LG(m,i,j)
    % LG  Local-to-Global index function.  Compute global
    % index k for unknown U_k, from local indices of (x_i,y_j).
    k = m * j + i;
    end % LG
