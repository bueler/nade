% TESTHEAT2D  Set up test case problem where we know exact solution,
% and call HEAT2D to solve numerically on sequence of refining grids.
% Measure and report ||E^h||_2 convergence rate.  Plot finest solution.

% test case definition:
uexact = @(x,y) x .* (1-x) .* cos((pi/2)*y);
f = @(x,y) - (2 + (pi^2/4) * x .* (1-x)) .* cos((pi/2)*y);

% measure ||E^h||_2 on refining grids
mm = [3 7 15 31 63 127];  % for h = 1/4, 1/8, 1/16, ...
%prepare to wait for finest:  mm = [3 7 15 31 63 127 255];
for s = 1:length(mm)
    h(s) = 1 / (mm(s) + 1);
    fprintf('solving case h=1/%d ...',mm(s)+1)
    tic
    [xx,yy,UU] = heat2d(mm(s),f);
    fprintf(' %.2f seconds\n',toc)
    UUex = uexact(xx,yy);
    err(s) = h(s) * norm(UU - UUex, 2);  % ||.||_2;  see page 252
end

% plot error analysis
figure(1)
p = polyfit(log(h),log(err),1);
errmodel = exp(p(2) + p(1)*log(h));
loglog(h,err,'o',h,errmodel,'r--')
axis tight
xlabel h,  ylabel('|E^h|_2')
title(sprintf('numerical error is |E^h|_2=O(h^{%.3f})',p(1)))

% plot final result
figure(2)
subplot(2,1,1)
surf(xx,yy,UU)
xlabel x, ylabel y, zlabel('numerical solution U_ij')
shading('flat')
subplot(2,1,2)
mesh(xx,yy,UU-UUex)
xlabel x, ylabel y, zlabel('error U_ij-u(x_i,y_j)')
