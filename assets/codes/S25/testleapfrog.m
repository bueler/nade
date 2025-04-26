% TESTLEAPFROG  Run convergence test on LEAPFROG.

a = -0.5;
tf = 10.0;
u0 = @(x) sin(8*pi*x);    % initial, and exact soln at tf

mm = [9 19 49 99 199 499];
hh = 1.0 ./ (mm + 1);     % = [0.1 0.05 0.02 0.01 0.005 0.002]
kk = hh;                  % refinement path

err = zeros(size(hh));
for s = 1:length(hh)
    x = (0:hh(s):1.0-hh(s))';
    if s == 1             % plot coarsest-case initial condition
        figure(1)
        xfine = linspace(0.0,1.0,501);
        plot(x,u0(x),'ko:',xfine,u0(xfine),'k')
        xlabel x
        legend('h=0.1 gridded initial cond','u(0,x) = sin(8 pi x)')
    end
    [x,U] = leapfrog(mm(s),a,kk(s),tf,u0(x));
    err(s) = norm(U - u0(x),'inf');
end

offset = 4;
figure(2)
loglog(hh,err,'ko')
hhp = hh(offset:end);  errp = err(offset:end);
p = polyfit(log(hhp),log(errp),1)
hold on,  loglog(hhp,exp(p(2) + p(1)*log(hhp)),'k--'),  hold off
text(0.003,1.1*err(4),['O(h^{' num2str(p(1)) '})'],'fontsize',16)
xlabel h,  ylabel('numerical error')
axis tight
set(gca(),'xtick',hh,...
          'xticklabel',{'0.1','0.05','0.02','0.01','0.005','0.002'})
grid on
