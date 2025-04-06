>> m = 25; h = 1/(m+1);
>> A = spdiags([ones(m,1), -2*ones(m,1), ones(m,1)], [-1 0 1], m, m);
>> A = (1/h^2) * A;
>> f = @(t,u) A * u;
>> [tt, uu] = ode45(f, [0, 10], ones(m,1));
