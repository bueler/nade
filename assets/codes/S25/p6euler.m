function yfinal = p6euler(m)
% P6EULER  Solves an ODE IVP on Assignment #1.  Does m steps.

% ODE IVP
f = @(t,Y) [Y(2);
            5 * Y(1) - 4 * Y(2)];
t0 = 2;  Y0 = [0; -1];

% Euler's method
tf = 4;  h = (tf - t0) / m;
t = t0;  Y = Y0;
for k = 1:m
    Y = Y + h * f(t,Y);
    t = t + h;
end
yfinal = Y(1);  % approximates y(tf) = y(4)
