% TAYLOR2 Compare Taylor O(k^2) method to Euler O(k^1).
% Also show Matlab's black-box (ode45) solution as exact.
% Solves:  u' = u + cos(u), u(0)=2.

f = @(t,u) u + cos(u);

% Euler solutions
figure(1),  hold on
for k = [1.0 0.5]
    t = 0.0:k:3.0;
    U = zeros(1,length(t));
    U(1) = 2.0;
    for j = 1:length(t)-1
        U(j+1) = U(j) + k * f(t(j),U(j));
    end
    plot(t,U,'bo:')
end

% Taylor O(k^2) solutions
df = @(t,u) 1.0 - sin(u);  % not needed by Euler, ode45
for k = [1.0 0.5]
    t = 0.0:k:3.0;
    U = zeros(1,length(t));
    U(1) = 2.0;
    for j = 1:length(t)-1
        fj = f(t(j),U(j));
        U(j+1) = U(j) + k * f(t(j),U(j)) * (1.0 + (k/2) * df(t(j),U(j)));
    end
    plot(t,U,'ro:')
end

% generate high-accuracy black box solution (and regard
% it as exact)
opt = odeset('RelTol',1.0e-12);
[tt,UU] = ode45(f,[0,3],2,opt);
plot(tt,UU,'k'),  hold off,  xlabel t
legend('Euler k=1','Euler k=0.5','Taylor2 k=1','Taylor2 k=0.5',...
       'exact (ode45)','Location','NorthWest','FontSize',12.0)
print -dpng taylor2.png
