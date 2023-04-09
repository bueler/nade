% ODERATES  Use test case x''+x=0, x(0)=1, x'(0)=0 to
% demonstrate expected convergence rates for final-time (tf=2)
% errors from explicit and non-adaptive solvers FEULER and RK4.

% problem
f = @(t,u) [u(2); -u(1)];
t0 = 0.0;
u0 = [1.0; 0.0];
tf = 2.0;

% measure final-time error in x(t) only
xexact = cos(tf);
erre = zeros(1,9);  errrk = erre;
N = 2.^(1:9);
for j = 1:9
    [te,ue] = feuler(f,u0,t0,tf,N(j));
    erre(j) = abs(ue(1,end) - xexact);
    [trk,urk] = rk4(f,u0,t0,tf,N(j));
    errrk(j) = abs(urk(1,end) - xexact);
end

% show measured error and fit convergence rates
h = 2.0 ./ N;
hfine = h(4:end);
loglog(h,erre,'bo','MarkerSize',10),  hold on
loglog(h,errrk,'ro','MarkerSize',10)
pe = polyfit(log(hfine),log(erre(4:end)),1);
loglog(hfine,exp(pe(1)*log(hfine)+pe(2)),'b:')
prk = polyfit(log(hfine),log(errrk(4:end)),1);
loglog(hfine,exp(prk(1)*log(hfine)+prk(2)),'r:')
hold off
xlabel h,  ylabel('final time error in x(t)')
legend(sprintf('Euler converges at O(k^{%.2f})',pe(1)),...
       sprintf('RK4 converges at O(k^{%.2f})',prk(1)),
       'FontSize',12.0,'Location','SouthEast')
axis tight
