% DCOMPARE  Compare FD formulas D_+, D_0 on variable frequency sines.

f = @(x,alf) sin(alf * x);
df = @(x,alf) alf .* cos(alf * x);

x = 1;
h = 0.001;
al = [1, 10, 100, 1000, 10000]'; 

Dperr = abs(df(x,al) - (f(x+h,al) - f(x,al))/h);
D0err = abs(df(x,al) - (f(x+h,al) - f(x-h,al))/(2*h));

format short e
[al, Dperr, D0err]
