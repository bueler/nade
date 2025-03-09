% RANDREALEIGS  Do a study of 3x3 random matrices to see
% what fraction have only real eigenvalues.

N = 100000;

countreal = 0;
countsymm = 0;
countsing = 0;
for m = 1:N
    A = randn(3,3);
    countreal = countreal + all(isreal(eig(A)));
    countsymm = countsymm + (norm(A - A') < 1.0e-10);
    countsing = countsing + (abs(det(A)) < 1.0e-10);
end

format short g
100 * countreal / N
100 * countsymm / N
100 * countsing / N
