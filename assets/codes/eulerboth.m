% EULERBOTH  Plot the solutions from applying (forward) Euler
% and backwards Euler methods to the simply ODE IVP
%   u'(t) = - 4 u(t),  u(0) = 2
% Step size k=1 is shown in Figure 1 and sizes k=0.5, 0.2 in
% Figure 2.

T = 4;

% exact soluton on fine grid
tt = linspace(0,T,501);
uexact = 2 * exp(-4 * tt);

% k = 1 step size
figure(1)
k = 1;
t = 0:k:T;
N = T / k;
uEULER = 2 * (-3).^(0:N);
uBE    = 2 * (1/5).^(0:N);
plot(t,uEULER,'o-',t,uBE,'o-')
hold on
plot(tt,uexact,'k','linewidth',1.0)
hold off
legend('Euler','backward Euler','exact solution')
xlabel t, ylabel('u(t)')
axis([0 T -7 8])
title('time step k=1')
print -dpng eulerbothone.png

% k = 0.5, 0.2 step sizes
figure(2)
hold on
for k = [0.5 0.2]
    t = 0:k:T;
    N = T / k;
    uEULER = 2 * (1 - 4*k).^(0:N);
    uBE    = 2 * (1/(1 + 4*k)).^(0:N);
    plot(t,uEULER,'bo-',t,uBE,'ro-')
end
hold on
plot(tt,uexact,'k','linewidth',2.0)
hold off
legend('Euler','backward Euler','Euler','backward Euler','exact solution')
xlabel t, ylabel('u(t)')
axis([0 2 -3 3])
title('time step sizes k=0.5, 0.2')
print -dpng eulerbothsmaller.png


% k = 0.1, 0.05 step sizes
figure(3)
hold on
for k = [0.1 0.05]
    t = 0:k:T;
    N = T / k;
    uEULER = 2 * (1 - 4*k).^(0:N);
    uBE    = 2 * (1/(1 + 4*k)).^(0:N);
    plot(t,uEULER,'bo-',t,uBE,'ro-')
end
hold on
plot(tt,uexact,'k','linewidth',2.0)
hold off
legend('Euler','backward Euler','Euler','backward Euler','exact solution')
xlabel t, ylabel('u(t)')
axis([0 2 -1 3])
title('time step sizes k=0.1, 0.05')
print -dpng eulerbothgood.png
