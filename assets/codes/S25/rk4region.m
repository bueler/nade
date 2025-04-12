% RK4REGION  Use Matlab's imagesc to plot the stability region of
% the classic RK4 scheme.

x = -3.5:0.01:3.5;
[xx, yy] = meshgrid(x, x);
zz = xx + i * yy;
RR = 1.0 + zz + 0.5 * zz.^2 + (1.0/6.0) * zz.^3 + (1.0/24.0) * zz.^4;
uu = (RR .* conj(RR) <= 1);
imagesc(x, x, 1 - uu, [-1, 1]),  colormap('gray')
grid on,  axis equal,  axis tight
