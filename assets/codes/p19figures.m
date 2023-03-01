% P19FIGURES  Reproduces Figure 2.4(b) and 2.5 in LeVeque (2007).

T = 2 * pi;  alpha = 0.7;  beta = 0.7;
m = 100;  h = T / (m+1);

figure(1),  clf
tt = h:h:T-h;
U0 = alpha + ((beta - alpha) / T) * tt;
pendulum(m,5,T,alpha,beta,U0);
title('Figure 2.4(b)')
print -dpng p19figure1.png

figure(2),  clf
U0 = 0.7 + sin(tt/2);
pendulum(m,8,T,alpha,beta,U0);
title('Figure 2.5')
print -dpng p19figure2.png
