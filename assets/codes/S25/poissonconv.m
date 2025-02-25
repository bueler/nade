% POISSONCONV  Show convergence of POISSON

% method of manufactured solutions
uexact = @(x,y) (x.^4 - x) .* (y.^4 - y);
f = @(x,y) 12.0 .* (x.^2 .* (y.^4 - y) + (x.^4 - x) .* y.^2);

% run and get error data
levs = 8;              % finest level is slow!
m = 2.^(1:levs) - 1;   % = [1, 3, 7, 15, ...]
h = 1.0 ./ (m+1);      % = [0.5, 0.25, 0.125, ...]
err = zeros(size(m));
for s = 1:levs
    fprintf('solving case h=1/%d ...',m(s)+1)
    tic
    [xx, yy, UU] = poisson(m(s), f);
    fprintf(' %.2f seconds\n', toc)
    Uex = uexact(xx, yy);
    err(s) = h(s) * norm(UU(:) - Uex(:), 2);
end

% generate figure with O(h^p) measured rate
loglog(h, err, 'ko')
p = polyfit(log(h), log(err), 1);
hold on,  loglog(h, exp(p(2) + p(1)*log(h)), 'k--'),  hold off
ostr = sprintf('O(h^{%.3f})', p(1))
text(1.5 * h(3), 1.5 * err(3), ostr, 'fontsize', 16)
xlabel('h', 'fontsize', 16)
axis tight
