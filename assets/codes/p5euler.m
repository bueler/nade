function Y = p5euler(m)
% P5EULER  Solves an ODE IVP on Assignment #1.  Does m steps.

f = @(t,Y) [Y(2);
            6*Y(1) - Y(2)];
t0 = 2;  Y0 = [0; -1];

tf = 4;
h = (tf - t0) / m;
t = t0;
Y = Y0;
for k = 1:m
    Y = Y + h * f(t,Y);
    t = t + h;
end
