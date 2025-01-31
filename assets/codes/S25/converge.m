% CONVERGE Calls SECOND and measures convergence rate in
% |.|_inf norm

f = @(x) exp(x);
mm = [5, 10, 20, 40, 100, 200, 400, 1000, 2000];
hh = 1 ./ (mm + 1);
for k = 1:length(mm)
    m = mm(k);
    [x,U] = second(f, 1, exp(1), m);
    err(k) = norm(U - exp(x), 'inf');
end
p = polyfit(log(hh), log(err), 1);
model = exp(p(2)) * exp(p(1) * log(hh));
loglog(hh, err, 'o', hh, model, 'r--')
axis tight
xlabel h,  ylabel('numerical error')
tag = sprintf('O(h^{%.3f})', p(1));
text(0.01, 1.0e-6, tag, 'color', 'r', 'fontsize', 24)
