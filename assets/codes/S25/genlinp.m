% GENLINP  Solve a 2-point boundary value problem
%   u''(x) + p u'(x) = 0,   u(0)=1,  u(1)=0
% numerically using centered FD, for an
% advection-dominating p value.  Visually compare
% to the exact solution.

p = -30;
% solve numerically and plot immediately
for m = [3 5 10 20 50 200 1000]
    fprintf('solving with m=%d ...',m)
    h = 1 / (m+1);  x = 0:h:1;
    A = -2 * speye(m);
    for j = 1:m
        if j > 1
            A(j,j-1) = -p * h / 2 + 1;
        end
        if j < m
            A(j,j+1) = p * h / 2 + 1;
        end
    end
    % for debugging:  full(A)
    A = (1/h^2) * A;
    F = zeros(m,1);
    F(1) = (1 / h^2) * (p * h / 2 - 1);
    U = A \ F;       % at interior points
    U = [1 U' 0];
    fprintf('done\n')
    labelstr = sprintf('.:;m=%d;',m);
    plot(x,U,labelstr,'MarkerSize',18.0),  hold on
end

% plot exact solution on a fine grid
xfine = 0:.001:1;
uexact = 1 - (1 - exp(-p*xfine)) ./ (1 - exp(-p));
plot(xfine,uexact,'k;u_{exact}(x);'),  hold off
xlabel x,  ylabel('u(x)')
legend('location','SouthWest')
axis tight
print -dpng genlinp.png
