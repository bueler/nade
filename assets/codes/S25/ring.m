% RING  Simulate the temperature in a ring of metal, a
% distributed and continuum object.  The model is the heat equation
% with additional advection, a partial differential equation (PDE)
% initial and boundary value problem:
%   u_t + v u_x = K u_xx + psi
% Here the heating source term psi(x) is only "on" over 0 < x < 1.
% The numerical solution is basic, but it is a reasonable first
% try for a finite difference (FD) method.  It combines first-order
% forward Euler in time with centered difference for the conduction
% term and simple upwinding for the advection.

u0 = 65.0;  % initial temperature; uniform
K = 0.3;    % diffusion (conduction over density * capacity) coefficient
v = 1.0;    % constant velocity; rotation
xc = 10.0;  % circumference of ring
psi = @(x,t) (x > 0.0 && x < 1.0) * 10.0;  % source term

tf = 10.0;   % final time; total duration
NN = 400;    % number of time steps
MM = 60;     % number of spatial steps

% set up for method
dt = tf / NN;  % time step
dx = xc / MM;  % spatial resolution
t = 0.0:dt:tf;
x = 0.0:dx:xc-dx;
R = v * dt / dx    % advection ratio
Q = K * dt / dx^2  % diffusion ratio

% fill source psi(x)
src = zeros(1, MM);
for j = 1:MM
    src(j) = psi(x(j), 0.0);
end

% FD method, with animation
u = zeros(1, MM);
u(:) = u0;
h = plot(x, u, 'o');
axis([0, xc, u0, u0 + 12.0])
ht = text(0.8 * xc, u0 + 10.0, '', 'fontsize', 20.0);
xlabel x,  ylabel('u(x,t) temperature'),  grid on
for n = 1:NN
    % plot current state
    pause(0.05)
    set(h, 'YData', u)
    set(ht, 'string', sprintf('t = %5.2f', t(n)))
    % the update
    ult = [u(MM), u(1:MM-1)];  % to the left, from periodicity
    urt = [u(2:MM), u(1)];     % to the right
    u = u - R * (u - ult) + Q * (urt - 2 * u + ult) + dt * src;
end
