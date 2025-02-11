function [h, Enorm] = verifygen(vcase, m)
% VERIFYGEN  Run a verification case on GENERALBVP.

uexact = @(x) sin(pi * x);
zerof = @(x) zeros(size(x));
if vcase == 1
    p = @(x) x.^2;
    q = zerof;
    f = @(x) - pi^2 * sin(pi * x) + pi * x.^2 .* cos(pi * x);
elseif vcase == 2
    p = zerof;
    q = @(x) x.^2;
    f = @(x) (x.^2 - pi^2) .* sin(pi * x);
else,  error('unknown vcase'),  end

h = (1 - 0) / (m + 1);
[x, U] = generalbvp(p, q, f, 0, 1, 0, 0, m);
Enorm = sqrt(h) * norm(U - uexact(x));
