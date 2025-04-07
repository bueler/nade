% STIFFCOMPARE  On linear systems u' = A u, where A is
% a finite difference approximation of the second
% derivative, so that the system is stiff, run two
% non-stiff Matlab ODE solvers (ode45, ode23) and
% compare a stiff solver (ode15s).

T = 2.0;
names = {'ode45 ', 'ode23 ', 'ode15s'};
fprintf('                k         N      run time (s)\n');
for m = [25, 50]
    h = 1/(m+1);
    v = ones(m,1);
    A = spdiags([v, -2*v, v], [-1, 0, 1], m, m);
    A = (1/h^2) * A;
    f = @(t,u) A * u;
    for method = 1:3
        tic
        if method == 1
            [tt, uu] = ode45(f, [0, T], v);
        elseif method == 2
            [tt, uu] = ode23(f, [0, T], v);
        else
            [tt, uu] = ode15s(f, [0, T], v);
        end
        runtime = toc;
        N = length(tt) - 1;
        kmax = max(diff(tt));
        fprintf('%s m=%3d:   %.1e %6d   %.2f\n', ...
                names{method}, m, kmax, N, runtime);
    end
end
