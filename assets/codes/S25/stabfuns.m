% STABFUNS  Plot stability functions R(x) relative to e^x.
% Note that R(z), for z complex, denotes the stability function
% of an ODE IVP scheme.  Here we just look along the real axis.

L = 4;
x = -L:.01:L;
Rfe = 1 + x;
Rem = 1 + x + (1/2) * x.^2;
Rrk = 1 + x + (1/2) * x.^2 + (1/6) * x.^3 + (1/24) * x.^4;
Rtr = (1 + (1/2) * x) ./ (1 - (1/2) * x);

plot(x, Rfe, x, Rem, x, Rrk, x, Rtr, x, exp(x), 'k')
xlabel('x = Re(z)'),  ylabel('R(x)')
legend('forward Euler', 'explicit midpoint', 'RK4', 'trapezoid', 'e^z', ...
       'location', 'south')
axis([-L, L, -5, 9]),  grid on
