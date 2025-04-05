% REGIONS  Plot stability regions for three second-order ODE IVP schemes
% This version uses imagesc instead of contour.

% set up
N = 300;
x = linspace(-5, 5, N+1);  y = x;
[xx,yy] = meshgrid(x, y);
zz = xx + i * yy;

% trapezoid
uu = ((abs(1 + zz/2) ./ abs(1 - zz/2)) <= 1);
subplot(1,3,1),  imagesc(x, y, 1 - uu, [-1, 1]),  colormap('gray')
grid on,  axis square,  title('TR')

% EM
uu = (abs(1 + zz + zz.^2/2) <= 1);
subplot(1,3,2),  imagesc(x, y, 1 - uu, [-1, 1]),  colormap('gray')
grid on,  axis square,  title('EM')

% AB2, both by pointwise decision and boundary locus
% note different scale
N = 200;
x = linspace(-2, 2, N+1);  y = x;
[xx,yy] = meshgrid(x, y);
zz = xx + i * yy;
uu = zeros(size(zz));
for j = 1:N+1
    for k = 1:N+1
        z = zz(j, k);
        zetas = roots([1, -(1+1.5*z), 0.5*z]);
        if abs(zetas(1)) < 1 && abs(zetas(2)) < 1
            uu(j, k) = 1.0;
        end
    end
end
subplot(1,3,3),  imagesc(x, y, 1 - uu, [-1, 1]),  colormap('gray')
th = 0:pi/100:2*pi;
zeta = exp(i * th);
zbl = (zeta.^2 - zeta) ./ (1.5 * zeta - 0.5);
hold on,  plot(real(zbl), imag(zbl), 'k', 'linewidth', 1.0),  hold off
grid on,  axis square,  title('AB2 (different scale)')
