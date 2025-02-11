% P13A  Do a verification study of GENERALBVP using the exact
% solution found in P13 (a) on Assignment 3.

zerof = @(x) zeros(size(x));
q = @(x) ones(size(x));
mm = [5, 10, 20, 40, 100, 200, 400, 1000, 2000];

for j = 1:length(mm)
    [x, U] = generalbvp(zerof, q, zerof, 0, 1, 2, 3, mm(j));
    uexact = 2 * cos(x) + (3 - 2 * cos(1)) * sin(x) / sin(1);
    h(j) = 1 / (mm(j) + 1);
    err(j) = norm(U - uexact, 'inf');
end
polyfit(log(h), log(err), 1)
