% TRBDF2REGION  Use Matlab's imagesc to plot the stability region of
% the TR-BDF2 scheme.  Compare RK4REGION.

figure(1)
x = -8:0.01:13;
y = -8:0.01:8;
[xx, yy] = meshgrid(x, y);
zz = xx + i * yy;
RR = (12 + 5 * zz) ./ ((3 - zz) .* (4 - zz));
uu = (abs(RR) <= 1);
imagesc(x, y, 1 - uu, [-1, 1]),  colormap('gray')
grid on,  axis tight

figure(2)
rrx = (12 + 5 * x) ./ ((3 - x) .* (4 - x));
plot(x, rrx, x, exp(x))
axis([-3 3 -2 3]),  grid on
xlabel x,  legend('R(x)', 'e^x')
