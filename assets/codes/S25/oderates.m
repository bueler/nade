% ODERATES  Use test case x''+4x=0, x(0)=1, x'(0)=0 to
% demonstrate expected convergence rates for final-time (tf=5)
% errors from explicit and non-adaptive solvers FE1 and RK4.

% problem
f = @(t,u) [u(2); -4*u(1)];
t0 = 0.0;  tf = 5.0;
u0 = [1.0; 0.0];

% measure final-time error in x(t) only
xexact = cos(2*tf);

% compare 8 levels of refinement
levs = 8;
errfe = zeros(1,levs);
errrk = errfe;
N = 10 * 2.^(1:levs);  % 20, 40, 80, ...
for j = 1:levs
    [tfe, ufe] = fe1(f,u0,t0,tf,N(j));
    errfe(j) = abs(ufe(1,end) - xexact);
    [trk, urk] = rk4(f,u0,t0,tf,N(j));
    errrk(j) = abs(urk(1,end) - xexact);
end

% show measured error and fit convergence rates
h = (tf - t0) ./ N;
loglog(h,errfe,'bo','MarkerSize',10),  hold on
loglog(h,errrk,'ro','MarkerSize',10)
pfe = polyfit(log(h),log(errfe),1);
loglog(h, exp(pfe(1)*log(h)+pfe(2)), 'b:')
prk = polyfit(log(h),log(errrk),1);
loglog(h, exp(prk(1)*log(h)+prk(2)), 'r:')
hold off,  xlabel h,  ylabel('final time error in x(t)')
legend(sprintf('Euler converges at O(k^{%.2f})',pfe(1)),...
       sprintf('RK4 converges at O(k^{%.2f})',prk(1)),
       'FontSize',12.0,'Location','SouthEast')
axis tight
