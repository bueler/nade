function heatcompare(m,N,tf,BEONLY)
% HEATCOMPARE Compute the forward and backward Euler (FE,BE)
% solutions of a heat equation problem
%   u_t = u_xx,  u(t,0) = u(t,1) = 0,  u(0,x) = f(x)
% Method-of-lines yields an linear ODE system
%   U'(t) = A U(t)
% Here U(t) is a column vector of size m.  We compute and report
% the maximum stable time step of FE, and estimated costs.
% Usage:  heatcompare(m,N,tf)
% where
%   m = number of space points
%   N = number of time steps
%   tf = final time.
% Example:
%   >> heatcompare(99,100,1.0)         % compute both FE and BE
%   >> heatcompare(999,1000,1.0,true)  % BE only at high resolution

if nargin < 4,  BEONLY = false;  end

% generate A and its eigenvalues
h = 1.0 / (m+1);
A = (1/h^2) * spdiags([ones(m,1), -2*ones(m,1), ones(m,1)],...
                      [-1, 0, 1], m, m);

% FE time restriction:  z = k lambda >= -2  for all eigenvalues
rhoA = max(abs(eig(A)));  % not all lambda are negative and real
maxkFE = 2.0 / rhoA;
kFE = 0.5 * maxkFE;       % safety factor of 0.5;  exact maxkFE
                          % does not generate actual decay
NFE = ceil(tf / kFE);

% initial condition  u(0,x) = f(x)
x = h:h:1-h;
U0 = zeros(m,1);
U0((x > 0.25) & (x < 0.5)) = 1.0;
subplot(3,1,1),  plot(x,U0,'k.','markersize',10)
grid on,  title('initial')

% compute FE
if ~BEONLY
    UFE = U0;
    for nn = 1:NFE
        UFE = UFE + kFE * A * UFE;  % FE step
    end
end

% compute BE; linear system is   (I - k A) U^{n+1} = U^n
kBE = tf / N;
UBE = U0;
for nn = 1:N
    UBE = (speye(m,m) - kBE * A) \ UBE;  % BE step
end

% exact solution: first two terms of Fourier sine series
uexact = (sqrt(2)/pi) * exp(-pi^2*tf) * sin(pi*x) ...
         + (1/pi) * exp(-4*pi^2*tf) * sin(2*pi*x);
errBE = norm(UBE - uexact','inf');
subplot(3,1,2:3),  hold on
plot(x,uexact,'r')
plot(x,UBE,'k.','markersize',10)
if ~BEONLY
    plot(x,UFE,'b.','markersize',10)
    legend('exact','from BE','from FE')
else
    legend('exact','from BE')
end
hold off,  xlabel x,  grid on,  title('final')

% results
fprintf('       stable time step for FE:  k = %.3e\n',kFE)
fprintf('        number of steps for FE:  N = %d\n',NFE)
fprintf('   time step used for N=%3d BE:  k = %.3e\n',N,kBE)
fprintf('  cost of FE (multiplications):  %d\n',NFE * 3 * m)
fprintf('  cost of BE (multiplications):  %d\n',N * 5 * m)
fprintf(' numerical error BE (inf norm):  %.3e\n',errBE)
