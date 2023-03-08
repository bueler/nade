function stepsizes(tol)
% STEPSIZES  Compare step sizes from ODE23 and ODE45 solutions
% of ODE IVP  u'(t) = 2 u(t) + exp(-t),  u(0) = 1
% Compare actual error at t=2 to requested RelTol.

if nargin < 1,  tol = 1.0e-4;  end

% problem
f = @(t,u) 2*u + exp(-t);
u0 = 1.0;  tf = 2.0;

% solve with black boxes and plot step sizes
[t23,u23] = ode23(f,[0,tf],u0,odeset('RelTol',tol,'AbsTol',tol));
[t45,u45] = ode45(f,[0,tf],u0,odeset('RelTol',tol,'AbsTol',tol));
plot(t23(2:end),diff(t23),'o',t45(2:end),diff(t45),'o')
xlabel t,  ylabel('adaptive time step k_n')
legend('ode23','ode45')
axis([0 tf 0 1.3*max(abs(diff(t45)))])

% compute U^N - uexact(tf)
uexact = @(t) (4/3) * exp(2*t) - (1/3) * exp(-t);
printf('for tol=%.3e:\n',tol)
printf('  ode23 takes %4d steps with final error |U^N-u(tf)|/|u(tf)| = %.3e\n',...
       length(t23)-1, abs(u23(end) - uexact(tf)) / abs(uexact(tf)))
printf('  ode45 takes %4d steps with final error |U^N-u(tf)|/|u(tf)| = %.3e\n',...
       length(t45)-1, abs(u45(end) - uexact(tf)) / abs(uexact(tf)))
