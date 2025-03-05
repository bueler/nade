% EXPLICIT12  Compare first- and second-order explicit ODE schemes
% on a long-time integration which shows their very different
% quality.  Forward Euler O(k^1) versus explicit midpoint O(k^2).
% ODE system is
%   v'(t) = w(t)
%   w'(t) = -9 v(t)
% with initial condition so that v(t) = cos(3 t) is exact solution.

v0 = 1.0;
w0 = 0.0;
f = @(u) [u(2); -9*u(1)];
u0 = [v0; w0];
uexact = @(t) cos(3 * t);

k = 0.02;  % time step
N = 1000;  % number of steps
t = 0:k:(N * k);
ufe = zeros(2, N+1);  % store whole soln
uem = zeros(2, N+1);  % store whole soln
ufe(:, 1) = u0;
uem(:, 1) = u0;
for j = 1:N
    ufe(:, j+1) = ufe(:, j) + k * f(ufe(:,j));
    ustar = uem(:, j) + (k / 2) * f(uem(:,j));
    uem(:, j+1) = uem(:, j) + k * f(ustar);
end

fprintf('result: ||u_FE - u_exact|| = %.2f\n', norm(ufe(1,:) - uexact(t)))
fprintf('        ||u_EM - u_exact|| = %.2f\n', norm(uem(1,:) - uexact(t)))

plot(t, ufe(1,:), t, uem(1,:), t, uexact(t))
xlabel t
ylabel('v(t)')
legend('forward Euler O(k^1)', 'explicit midpoint O(k^2)', 'exact', ...
       'location', 'northwest')
