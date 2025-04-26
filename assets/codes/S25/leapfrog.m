function [x, U] = leapfrog(m, a, k, tf, u0)
% LEAPFROG  Use the centered-space and midpoint-time FD method
% on the advection equation
%     u_t + a u_x = 0
% for u(t,x), with periodic boundary conditions on [0,1] and
% initial condition u(x,0) = u0(x).  Uses m+1 grid points and a
% spatial grid with h = 1/(m+1).  Forms the matrix A in the MOL system
%     U'(t) = A U(t).
% Uses time steps k and final time tf.  First step is by
% RK2 which is also second order.
% Usage:  [x, U] = leapfrog(m, a, k, tf, u0)
% See TESTLEAPFROG.

% generate space grid and time-step
h = 1.0 / (m+1);
x = (0:h:1.0-h)';              % interpret as periodic grid
NN = ceil(tf / k);             % satisfy  k NN = tf
k = tf / NN;

% sparse circulant matrix A
A = spdiags([-ones(m+1,1), zeros(m+1,1), ones(m+1,1)],...
            [-1, 0, 1],m+1,m+1);
A(m+1,1) = 1.0;  A(1,m+1) = -1.0;
A = - (a/(2.0*h)) * A;
% full(A), spy(full(A))        % check A (for m small)

% solve advection equation
Uold = u0(:);                  % user's initial condition as a column
U = Uold + (k/2) * A * Uold;   % first step is RK2, eqn (5.30)
U = Uold + k * A * U;
for n = 2:NN
    Unew = Uold + 2.0 * k * A * U; % leapfrog is multi-step
    Uold = U;                      %   so update two vectors
    U = Unew;
end
