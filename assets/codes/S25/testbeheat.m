function testbeheat
% TESTBEHEAT  Run a convergence test on BEHEAT.

D = 1/20;
tf = 0.1;
u0 = @(x) sin(5*pi*x);
uexact = @(x,t) exp(-25*pi^2*D*t) * sin(5*pi*x);
mm = [20-1, 100-1, 200-1, 500-1, 1000-1, 2000-1];
hh = 1.0 ./ (mm + 1);
kk = hh;                         % refinement path

err = zeros(size(hh));
for j = 1:length(mm)
    fprintf('case h=%.4f', hh(j))
    [x,U] = beheat(D, tf, u0, mm(j), kk(j));
    fprintf(' ... done\n')
    err(j) = norm(U - uexact(x, tf), 'inf');
end
loglog(hh, err, 'ko'),  hold on
p = polyfit(log(hh), log(err), 1);
loglog(hh, exp(p(2) + p(1)*log(hh)), 'k--')
text(0.003, err(4), ['O(h^{' num2str(p(1)) '})'], ...
        'fontsize', 14)
set(gca(), 'xtick', [0.05 0.01 0.005 0.001 0.0005], ...
           'xticklabel', {'0.05','0.01','0.005','0.001','0.0005'})
xlabel h,  ylabel('numerical error'),  axis tight,  grid on
