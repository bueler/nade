% EDGEOFSTABILITY  Show how a linear system behaves
% when forward Euler is used on either side of the
% maximum stable time step k_max.  Note
%   >> A = [1 -2 3; 4 -5 6; 0 1 -2]; eig(A)
% gives -5, -1, 0, and thus theory shows
%   kstab = 2/5 = 0.4

f = @(t,u) [u(1)-2*u(2)+3*u(3); 4*u(1)-5*u(2)+6*u(3); u(2)-2*u(3)];
eta = [4, 9, -3]';
tf = 10;

subplot(5,1,1:2)
[tt, uu] = fe1(f, eta, 0, tf, 40);
plot(tt, uu),  title('k=0.25')
[ttt, uuu] = fe1(f, eta, 0, tf, 400);  % near exact
hold on,  plot(ttt, uuu, 'k', 'linewidth', 1),  hold off
subplot(5,1,3)
[tt, uu] = fe1(f, eta, 0, tf, 30);
plot(tt, uu),  title('k=0.33')
subplot(5,1,4)
[tt, uu] = fe1(f, eta, 0, tf, 25);
plot(tt, uu),  title('k=0.40')
subplot(5,1,5)
[tt, uu] = fe1(f, eta, 0, tf, 20);
plot(tt, uu),  title('k=0.50')
xlabel t
