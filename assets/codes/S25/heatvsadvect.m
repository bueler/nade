function heatvsadvect(eqn, ic)
% HEATVSADVECT  Compare numerical methods for heat equation
% versus advection:
%   u_t = D u_xx     heat  (with D=1)
%   u_t + a u_x = 0  advection  (with a=1)
% Uses periodic boundary conditions on [0,1] for both problems,
% and it considers two cases for initial condition: pure wave,
% square.  Numerical method is FTCS for heat and first-order
% upwinding for advection.

h = 0.01;  % grid spacing
tf_heat = 0.02;
D_heat = 1.0;
tf_adv = 0.2;
a_adv = 1.0;

figure  % open new figure
MM = round(1.0 / h);
x = linspace(0.0, 1.0 - h, MM);  % as periodic grid
if ic == 1
    u = sin(4 * pi * x);
    bottom = -1.2;  % for plotting
else
    u = zeros(size(x));
    mark = (x > 0.3) & (x < 0.6);
    u(mark) = 1.0;
    bottom = -0.2;  % for plotting
end

if eqn == 1  % heat
    tf = tf_heat;
    k = 0.5 * h^2 / (2 * D_heat);  % stability criterion for FTCS; half of limit
    NN = ceil(tf / k);
    k = tf / NN;
    plot(x, u);
    axis([0, 1, bottom, 1.2]),  hold on,  grid on
    xlabel x
    r = D_heat * k / h^2;
    fprintf('simulating N=%d steps of heat equation on ic %d ...\n', NN, ic)
    for n = 1:NN
        uleft = [u(end), u(1:end-1)];
        uright = [u(2:end), u(1)];
        u = u + r * (uright - 2 * u + uleft);
    end
    plot(x, u);  hold off
    legend('u(0.000,x)', sprintf('u(%.3f,x)', tf))
    title(sprintf('heat equation, i.c. case %d', ic))
else  % advection
    tf = tf_adv;
    k = 0.5 * h / abs(a_adv);  % stability criterion for upwind; half of CFL limit
    NN = ceil(tf / k);
    k = tf / NN;
    plot(x, u);
    axis([0, 1, bottom, 1.2]),  hold on,  grid on
    xlabel x
    r = a_adv * k / h;
    fprintf('simulating N=%d steps of advection equation on ic %d ...\n', NN, ic)
    if a_adv >= 0
        for n = 1:NN
            uleft = [u(end), u(1:end-1)];
            u = u - r * (u - uleft);
        end
    else
        for n = 1:NN
            uright = [u(2:end), u(1)];
            u = u - r * (uright - u);
        end
    end
    plot(x, u);  hold off
    legend('u(0.000,x)', sprintf('u(%.3f,x)', tf))
    title(sprintf('advection equation, i.c. case %d', ic))
end
