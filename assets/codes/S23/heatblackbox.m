% HEATBLACKBOX  Demo the advantage of an MOL discretization
% by showing how easy it is to generate a very accurate
% solution to same problem as in P34 on Assignment #8.

m = 99;  h = 1/(m+1);
A = -2 * speye(m,m);
for j = 1:m-1,  A(j,j+1)=1.0;  A(j+1,j)=1.0;  end
D = 1/20;  A = (D/h^2) * A;
f = @(t,U) A*U;
x = h:h:1-h;  eta = sin(5*pi*x);
[tt,UU] = ode23(f,[0.0,0.1],eta');
plot(x,UU(end,:),'-o'),  xlabel x
hold on, plot(x,exp(-25*pi^2*D*0.1)*eta,'r'),  hold off

legend('U(tf)','exact solution u(tf,x)')
#print -dpng heatblackbox.png
