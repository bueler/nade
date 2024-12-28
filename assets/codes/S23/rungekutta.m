% RUNGEKUTTA Compare O(k^1), O(k^2), and O(k^4) Runge-Kutta
% methods, i.e. Euler, ET, RK4.
% Also show Matlab's black-box (ode45) solution as exact.
% Solves:  u' = u + cos(u), u(0)=2.

f = @(t,u) u + cos(u);
k = 1.0;
t = 0.0:k:3.0;

% Euler
U = zeros(1,length(t));
U(1) = 2.0;
for j = 1:length(t)-1
    U(j+1) = U(j) + k * f(t(j),U(j));
end
plot(t,U,'bo:'),  hold on

% explicit trapezoid (ET)
U = zeros(1,length(t));
U(1) = 2.0;
for j = 1:length(t)-1
    Ustar = U(j) + k * f(t(j),U(j));
    U(j+1) = U(j) + (k/2) * ( f(t(j),U(j)) + f(t(j+1),Ustar) );
end
plot(t,U,'ro:')

% RK4
U = zeros(1,length(t));
U(1) = 2.0;
for j = 1:length(t)-1
    Y1 = U(j);
    f1 = f(t(j),Y1);
    Y2 = U(j) + (k/2) * f1;
    f2 = f(t(j)+k/2,Y2);
    Y3 = U(j) + (k/2) * f2;
    f3 = f(t(j)+k/2,Y3);
    Y4 = U(j) + k * f3;
    f4 = f(t(j+1),Y4);
    U(j+1) = U(j) + (k/6) * ( f1 + 2 * f2 + 2 * f3 + f4 );
end
plot(t,U,'go:')

% generate high-accuracy black box solution (and regard
% it as exact)
opt = odeset('RelTol',1.0e-12);
[tt,UU] = ode45(f,[0,3],2,opt);
plot(tt,UU,'k'),  hold off,  xlabel t
legend('Euler','explicit trapezoid','RK4',...
       'exact (ode45)','Location','NorthWest','FontSize',12.0)
print -dpng rungekutta.png
