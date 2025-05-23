% FDRATES  Reproduce Figure 1.2 in LeVeque showing errors
% from several FD formulas.

u = @sin;                         % u is a function handle
x = 1.0;   exact = cos(x);
h = [0.1 0.05 0.01 0.005 0.001];
Dplus = (u(x+h) - u(x)) ./ h;
D0    = (u(x+h) - u(x-h)) ./ (2.0 * h);
D3    = (2.0 * u(x+h) + 3.0 * u(x) - 6.0 * u(x-h) + u(x-2.0*h)) ...
        ./ (6.0 * h);
loglog(h,abs(Dplus - exact),'k.-','markersize',18),  hold on
loglog(h,abs(D0 - exact),'k.-','markersize',18)
loglog(h,abs(D3 - exact),'k.-','markersize',18)
text(0.0015,2.0e-3,'D_+','fontsize',20)
text(0.0015,1.0e-6,'D_0','fontsize',20)
text(0.0015,1.0e-9,'D_3','fontsize',20)
hold off,  axis([5.0e-4 0.2 1.0e-11 0.1])
