function [x, U] = beheat(D,tf,u0,m,k)
% BEHEAT Use backward Euler and centered-space FD on the heat eqn
%     u_t = D u_xx
% for u(x,t), with boundary conditions u(0,t)=u(1,t)=0 and initial
% condition u(x,0)=u0(x).  Uses m interior points and spatial
% grid with h = 1/(m+1) to form the MOL matrix A.  Uses time steps
% k and final time tf to solve the problem.
% Usage:  [x, U] = beheat(D,tf,u0,m,k)
% Example: TESTHEATS.

% generate space grid, matrix A, and time-step
h = 1.0 / (m+1);
x = (h:h:1.0-h)';         % column vector
A = (D/h^2) * spdiags([ones(m,1), -2*ones(m,1), ones(m,1)],...
                      [-1, 0, 1],m,m);
NN = ceil(tf / k);        % assure:  k NN = tf
k = tf / NN;

U = u0(x);                % evaluate user's initial condition
B = speye(m,m) - k * A;   % matrix for taking time step
for n = 1:NN
    U = B \ U;            % the backward Euler step
end
