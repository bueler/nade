% TESTLINPEND  Generate convergence figure for LINPEND.

% test case
a = 0;  b = 1;  alpha = 2;  beta = 3;
uexact = @(x) ((3 - 2 * cos(1)) / sin(1)) * sin(x) + 2 * cos(x);
mm = [5 10 20 40 80 160 320 640];
for k = 1:length(mm)
    h(k) = (b - a) / (mm(k) + 1);
    [x,U] = linpend(mm(k),a,b,alpha,beta);
    err(k) = sqrt(h(k)) * norm(U - uexact(x));
end
p = polyfit(log(h),log(err),1);
loglog(h,err,'o',h,exp(p(2) + p(1)*log(h)),':')
title(sprintf('convergence at rate O(h^{%.3f})',p(1)),'fontsize',18)
xlabel('h','fontsize',14),  ylabel('|U - uexact|_2','fontsize',14)
axis tight
