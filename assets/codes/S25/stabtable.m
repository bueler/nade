% STABTABLE Nearly reproduce Table 7.1 in LeVeque.  Calls FE1.

lambda = -2100;
f = @(t,u) lambda * (u - cos(t)) - sin(t);
u0 = 1.0;
T = 2.0;

klist = [0.001, 0.000976, 0.000950, 0.000800, 0.000400];
% what exactly did LeVeque use?
% Nlist = [2000, 2050, 2106, 2500, 5000];
for j = 1:length(klist)
    k = klist(j);
    N = ceil(T / k);
    [tt, uu] = fe1(f, u0, 0.0, T, N);
    fprintf('%.6f:  %.5e\n', T/N, abs(uu(end) - cos(T)))
end
