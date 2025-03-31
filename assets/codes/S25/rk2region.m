% RK2REGION  Use Matlab's contourf to plot the stability
% region of the explicit midpoint rule, an RK2 scheme.

x = -3:0.1:3;
[xx, yy] = meshgrid(x, x);
rr = (1 + xx + 0.5*xx.^2 - 0.5*yy.^2).^2 + (yy + xx .* yy).^2;
contourf(xx, yy, rr, [0 1])
colormap(gray)
grid on,  axis equal,  axis tight
