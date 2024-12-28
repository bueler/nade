% VERIFYSH  Verify STEADYHEAT by measuring errors on refining grids.

f = @(x) 2 - 12 * x + 12 * x.^2;
uexact = @(x) x.^2 .* (1 - x).^2;
mlist = [10 20 40 80 160 320 640 1280];
hlist = 1 ./ (mlist + 1);
for k = 1:length(mlist)
    [x,U] = steadyheat(mlist(k),f,0.0,0.0);
    % see p. 17 in LeVeque about this norm
    err(k) = sqrt(hlist(k)) * norm(U - uexact(x),2);
end
p = polyfit(log(hlist),log(err),1);
fprintf('convergence at rate O(h^p) with p = %f\n',p(1))
loglog(hlist,err,'o',hlist,exp(p(2) + p(1)*log(hlist)),'r--')
xlabel h,  ylabel('numerical error')
text(0.01,0.0002,sprintf('O(h^{%.3f})',p(1)),...
     'Color','r','FontSize',14)
axis tight
