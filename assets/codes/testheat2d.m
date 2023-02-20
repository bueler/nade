function testheat2d(m)
% TESTHEAT2D  Set up test case problem where we know
% exact solution,
%   u(x,y) = x (1-x) cos((pi/2) y)
% and call HEAT2D to solve numerically.  Then measure
% numerical error.

f = @(x,y) - (2 + (pi^2/4) * x .* (1-x)) .* cos((pi/2)*y);
[xx,yy,UU] = heat2d(m,f);

uexact = @(x,y) x .* (1-x) .* cos((pi/2)*y);
UUex = uexact(xx,yy);
% arrays to debug:
%xx, yy, UU, UUex

norm(UU-UUex,'inf')
