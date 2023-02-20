function [xx,yy,UU] = heat2d(m)
% HEAT2D  Solve on the unit square [0,1]^2:
%   u_xx + u_yy = f(x,y)
% with u=0 on x=0, x=1, and y=0 sides, but u_y=0 on the y=1 side.
% Example:
%    >> [xx,yy,UU] = heat2d(2)

h = 1 / (m+1);
[xx,yy] = meshgrid(0:h:1,0:h:1);

n = m * (m+1);   % number of unknowns
A = zeros(n,n);
F = zeros(n,1);
for i = 1:m
    for j = 1:m+1
        k = LG(m,i,j)
        A(k,k) = -4;
        if i > 1
            A(k,LG(m,i-1,j)) = 1;
        end
        if i < m
            A(k,LG(m,i+1,j)) = 1;
        end
        if j > 1 && j < m + 1
            A(k,LG(m,i,j-1)) = 1;
        elseif j == m + 1
            A(k,LG(m,i,j-1)) = 2;
        end
        if j < m + 1
            A(k,LG(m,i,j+1)) = 1;
        end
    end
end
spy(A)
end % heat2d

    function k = LG(m,i,j)
    k = m * (j - 1) + i;
    end % LG
