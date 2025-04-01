% VERIFYBDF2  Use test case x''+4x=0, x(0)=1, x'(0)=0 to
% demonstrate expected convergence rate for final-time (tf=5)
% errors from the second-order implicit solver BDF2.

% problem
A = [0, 1; -4, 0];
g = @(t) zeros(2,1);
t0 = 0.0;  tf = 5.0;
u0 = [1.0; 0.0];
xexact = cos(2*tf);

% compute 9 levels of refinement
levs = 9;
err = zeros(1, levs);
N = 10 * 2.^(1:levs);  % 10, 20, 40, 80, ...
for j = 1:levs
    [tt, zz] = bdf2(A, g, u0, t0, tf, N(j));
    err(j) = abs(zz(1,end) - xexact);
end

% show measured error, but fit only the finer-mesh errors
% when computing the convergence rate
h = (tf - t0) ./ N;
loglog(h, err, 'bo', 'MarkerSize', 10),  hold on
p = polyfit(log(h(4:end)), log(err(4:end)), 1);
loglog(h(4:end), exp(p(1) * log(h(4:end)) + p(2)), 'r:')
hold off,  xlabel h,  axis tight
ylabel('final time error in x(t)')
legend(sprintf('BDF2 converges at O(k^{%.2f})', p(1)),
       'FontSize',12.0,'Location','SouthEast')
