function [t,x,T] = heat(m,N)
% HEAT  Demonstrate the explicit, first-order Euler, centered-in-space
% scheme for the heat equation  T_t = K T_xx.  This code fixes the
% space interval [0,L], the time interval [0,tau], the diffusivity K
% and the initial condition T(0,x) = g(x).  Boundary conditions are
% periodic.  User inputs set the number of subintervals in space(m)
% and time (N).  Produces graphable output suitable for surface plot.
% Warning 1: Choosing the "wrong" m and N will reveal instability!
%            (See heat(30,20).)  Course content will address this.
% Warning 2: There is no verification case here.  Why should we trust
%            the output?
% Warning 3: This code saves everything as arrays.  This is not a
%            good design for larger simulations.
% Example:
%   >> [t,x,T] = heat(20,20);
%   >> surf(x,t,T),  shading flat
%   >> xlabel('x'),  ylabel('t'),  zlabel('temperature T(t,x)')

% set parameters (Warning 4: No input checking!)
K = 0.01;
L = 1.0;
tau = 1.0;
h = L / m;
delta = tau / N;

% create mesh:  (t,x) in [0,tau] x [0,L]
[t,x] = ndgrid(0:delta:tau,h/2:h:L-h/2);

% initial condition:  g(x) = 1 on [0.1,0.3], g(x) = 0 otherwise
g = ones(size(x(1,:)));
g(x(1,:) < 0.1) = 0.0;
g(x(1,:) > 0.3) = 0.0;

% time-stepping
T = zeros(size(x));    % allocate space for solution
T(1,:) = g;            % set initial condition
R = K * delta / h^2;   % the key parameter for stability!
jj = 2:m-1;            % interior node indices; thus m >= 3 required
for n = 1:N            % step forward in time to fill T(t,x)
    T(n+1,1)  = T(n,1)  + R * ( T(n,2) - 2 * T(n,1) + T(n,m) );   % periodic
    T(n+1,jj) = T(n,jj) + R * ( T(n,jj+1) - 2 * T(n,jj) + T(n,jj-1) );
    T(n+1,m)  = T(n,m)  + R * ( T(n,1) - 2 * T(n,m) + T(n,m-1) ); % periodic
end
