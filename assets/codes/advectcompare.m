% ADVECTCOMPARE  Compare
%   Lax-Friedrichs
%   first-order upwind
%   leapfrog
%   Lax-Wendroff
% in single figure, on
%   u_t + 0.5 u_x = 0
% over t,x in [0,2]x[0,1] with square initial condition
% and periodic boundary conditions.
% See also: ADVECTTRY
% Also try:  tf = 6.0
%            square = false
%            m = 200

% problem
a = 0.5;
tf = 2.0;
square = true;
m = 80;

% space-time grid
h = 1.0 / m;
k = 0.9 * (h / abs(a));   % safety factor: nu = ak/h = 0.9 < 1
NN = ceil(tf / k);
k = tf / NN;
nu = a * k / h;           % Courant number
x = 0.0:h:1.0-h;

% initial condition and exact solution
if square
    u0 = ones(size(x));
    u0(x < 0.4) = 0.0;
    u0(x > 0.6) = 0.0;
else
    u0 = sin(3*pi*x).^2;
end
uexact = u0;              % tf must be even integer!

% Lax-Friedrichs:  (10.6) in LeVeque
figure(1)
U = u0;
for n = 1:NN
    Uright = [U(2:m) U(1)];
    Uleft  = [U(m) U(1:m-1)];
    Uav = 0.5 * (Uleft + Uright);
    U = Uav - (nu/2) * (Uright - Uleft);
    if n == 1
        U1lf = U;
    end
end
Ulf = U;
subplot(4,1,1)
plot(x,Ulf,'o-',x,uexact,'r',linewidth=0.5)
axis([0 1 -0.5 1.5])
xlabel x
legend('Lax-Friedrichs',sprintf('exact u(%.2f,x)',tf))
text(0.02,1.2,sprintf('|U-uexact|_1 = %.3f',h*norm(Ulf - uexact,1)))
title(sprintf('m=%d',m))

% upwind:  (10.21) in LeVeque
U = u0;
for n = 2:NN
    Uleft  = [U(m) U(1:m-1)];
    U = U - nu * (U - Uleft);   % assumes a >= 0
end
Uup = U;
subplot(4,1,2)
plot(x,Uup,'o-',x,uexact,'r',linewidth=0.5)
axis([0 1 -0.5 1.5])
xlabel x
legend('upwind',sprintf('exact u(%.2f,x)',tf))
text(0.02,1.2,sprintf('|U-uexact|_1 = %.3f',h*norm(Uup - uexact,1)))

% leapfrog:  (10.13) in LeVeque
% with first step from LF
Uold = u0;
U = U1lf;
for n = 2:NN
    Uright = [U(2:m) U(1)];
    Uleft  = [U(m) U(1:m-1)];
    Unew = Uold - nu * (Uright - Uleft);
    Uold = U;
    U = Unew;
end
Uleap = U;
subplot(4,1,3)
plot(x,Uleap,'o-',x,uexact,'r',linewidth=0.5)
axis([0 1 -0.5 1.5])
xlabel x
legend('leapfrog',sprintf('exact u(%.2f,x)',tf))
text(0.02,1.2,sprintf('|U-uexact|_1 = %.3f',h*norm(Uleap - uexact,1)))

% Lax-Wendroff: (10.18) in LeVeque
U = u0;
for n = 1:NN
    Uright = [U(2:m) U(1)];
    Uleft  = [U(m) U(1:m-1)];
    U = U - (nu/2) * (Uright - Uleft) + (nu^2/2) * (Uright - 2 * U + Uleft);
end
Ulw = U;
subplot(4,1,4)
plot(x,Ulw,'o-',x,uexact,'r',linewidth=0.5)
axis([0 1 -0.5 1.5])
xlabel x
legend('Lax-Wendroff',sprintf('exact u(%.2f,x)',tf))
text(0.02,1.2,sprintf('|U-uexact|_1 = %.3f',h*norm(Ulw - uexact,1)))
