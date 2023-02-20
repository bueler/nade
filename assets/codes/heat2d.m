function [xx,yy,UU] = heat2d(m,fsource)
% HEAT2D  Solve on the unit square [0,1]^2:
%   u_xx + u_yy = f(x,y)
% with Dirichlet conditions u=0 on the x=0, x=1, and y=1 sides,
% but u_y=0 on the y=0 side.  The grid has local (left) indices
% and global (right) indices as follows in the m = 3 case, with
% Dirichlet boundary indicated by zero:
%   j=4  0   0   0   0   0           0   0   0   0   0
%   j=3  0   *   *   *   0           0   10  11  12  0
%   j=2  0   *   *   *   0           0   7   8   9   0
%   j=1  0   *   *   *   0           0   4   5   6   0
%   j=0  0   *   *   *   0           0   1   2   3   0
%       i=0 i=1 i=2 i=3 i=4                  k
% Note (x_i,y_j) = (i h,j h), where h = /(m+1), is grid location.
% Usage:
%    >> [xx,yy,UU] = heat2d(m,fsource)
% Example:  See TESTHEAT2D.

h = 1 / (m+1);
[xx,yy] = ndgrid(0:h:1,0:h:1);

n = m * (m+1);             % number of unknowns
A = zeros(n,n);
F = zeros(n,1);
%VV = zeros(size(xx));     % for debugging
for j = 0:m
    for i = 1:m
        k = LG(m,i,j);     % at (x_i,y_j); will fill kth row of A
        %VV(i+1,j+1) = k;  % for debugging
        F(k) = fsource(xx(i+1,j+1), yy(i+1,j+1));
        A(k,k) = -4;
        if i > 1
            A(k,LG(m,i-1,j)) = 1;
        end
        if i < m
            A(k,LG(m,i+1,j)) = 1;
        end
        if j > 0
            A(k,LG(m,i,j-1)) = 1;
        end
        if j == 0
            A(k,LG(m,i,j+1)) = 2;
        elseif j < m
            A(k,LG(m,i,j+1)) = 1;
        end
    end
end
A
A = (1/h^2) * A;
% to check matrix structure:
%spy(A)

% for debugging:  to check grid ordering, uncomment
% all "VV" lines above and here:
%VV, xx, yy

% solution as an n x 1 vector
U = A \ F;

% fill solution into grid positions
UU = zeros(size(xx));
UU(2:m+1,1:m+1) = reshape(U,m,m+1);

end % heat2d

    function k = LG(m,i,j)
    k = m * j + i;
    end % LG
