function [x, U] = beheat(D, tf, u0, m, k)
% BEHEAT Use centered FD in space, and backward Euler in time, on
% the heat equation
%     u_t = D u_xx
% for u(x,t).  The boundary conditions are u(0,t)=u(1,t)=0, and
% the initial condition is u(x,0)=u0(x).  Uses m interior points
% and spatial a grid with h = 1/(m+1) to form the MOL matrix A.
% Uses time steps k and final time tf to solve the problem.
% Returns only the final-time state U ~ u(x,tf).
% Usage:  [x, U] = beheat(D, tf, u0, m, k)
% Example: TESTBEHEAT.

% generate space grid, matrix A, and time-step
h = 1.0 / (m+1);
x = (h:h:1.0-h)';         % column vector
A = (D/h^2) * spdiags([ones(m,1), -2*ones(m,1), ones(m,1)],...
                      [-1, 0, 1],m,m);
N = ceil(tf / k);        % force to be true:  k NN = tf
k = tf / N;

U = u0(x);                % evaluate user's initial condition
B = speye(m,m) - k * A;   % matrix for taking time step
for q = 1:N
    U = B \ U;            % the backward Euler step
end
