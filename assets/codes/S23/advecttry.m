% ADVECTTRY  Demonstrate FTCS and Lax-Friedrichs on
%   u_t + 0.5 u_x = 0
% over t,x in [0,2]x[0,1] with square initial condition
% and periodic boundary conditions.  Generates three
% figures.

a = 0.5;
tf = 2.0;
for m = [20 80 320]
    h = 1.0 / m;
    k = 0.9 * (h / abs(a));   % safety factor: nu = ak/h = 0.9 < 1
    % k = h / abs(a);         % suspicious result?
    NN = ceil(tf / k)
    k = tf / NN;
    nu = a * k / h;           % Courant number
    x = 0.0:h:1.0-h;
    u0 = ones(size(x));
    u0(x < 0.4) = 0.0;
    u0(x > 0.6) = 0.0;
    uexact = u0;              % why? (note tf matters!)
    % FTCS
    U = u0;
    for n = 1:NN
        Uright = [U(2:m) U(1)];
        Uleft  = [U(m) U(1:m-1)];
        U = U - (nu/2) * (Uright - Uleft);
    end
    Uftcs = U;
    % Lax-Friedrichs
    U = u0;
    for n = 1:NN
        Uright = [U(2:m) U(1)];
        Uleft  = [U(m) U(1:m-1)];
        Uav = 0.5 * (Uleft + Uright);
        U = Uav - (nu/2) * (Uright - Uleft);
    end
    Ulf = U;
    norm(Uftcs - uexact,2)
    norm(Ulf - uexact,2)
    figure()
    plot(x,Uftcs,'o-',x,Ulf,'o-',x,uexact,'r')
    axis([0 1 -0.5 1.5])
    xlabel x
    legend('FTCS','Lax-Friedrichs',sprintf('u(%.2f,x)',tf))
    title(sprintf('m=%d',m))
end
