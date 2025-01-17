% DERRORS  Plot errors in finite differences approximations
% of f'(x), in a particular case.

f = @(x) exp(x);
df = f;
x = 3.0;
h = 10.^(-(0:7));

Dp = @(f, x, h) (f(x+h) - f(x)) / h;
Dm = @(f, x, h) (f(x) - f(x-h)) / h;
D0 = @(f, x, h) (f(x+h) - f(x-h)) / (2*h);

for j = 1:8
    Dperr(j) = abs(Dp(f,x,h(j)) - df(x));
    Dmerr(j) = abs(Dm(f,x,h(j)) - df(x));
    D0err(j) = abs(D0(f,x,h(j)) - df(x));
end

loglog(h, Dperr, 'o-', h, Dmerr, 'o-', h, D0err, 'o-')
xlabel h,  ylabel('absolute error')
legend('D_+ f', 'D_- f', 'D_0 f', ...
       'location', 'northwest')
grid on
