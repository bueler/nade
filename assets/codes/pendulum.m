function pendulum(m,K,T,alpha,beta,U0)
% PENDULUM  Solve ODEBVP for nonlinear pendulum problem,
%     U''(t) = - sin(U(t)),
% where U(t) is the pendulum angle from the vertical.
% Uses boundary conditions U(0) = alpha, U(T) = beta.
% Sets-up and solves a nonlinear system of equations
%     G(U) = 0
% by K iterations of Newton's method, following section
% 2.16 in LeVeque (2007).  Plots and labels the iterates.
% Usage:
%     pendulum(m,K,T,alpha,beta,U0)
% where
%     m = number of interior points in grid
%     K = number of Newton iterations
%     T = length of time interval
%     alpha = U(0),  beta = U(T)
%     U0 = initial iterate
% Example: P19FIGURES

h = T / (m+1);

    function z = G(U)
        % Evaluate residual in nonlinear equations.
        z = zeros(m,1);
        z(1) = (alpha - 2.0 * U(1) + U(2)) / h^2 + sin(U(1));
        for j = 2:m-1
            z(j) = (U(j-1) - 2.0 * U(j) + U(j+1)) / h^2 + sin(U(j));
        end
        z(m) = (U(m-1) - 2.0 * U(m) + beta) / h^2 + sin(U(m));
    end

    function A = Jacobian(U)
        % Evaluate Jacobian in nonlinear equations.
        A = diag(-2.0 / h^2 + cos(U));
        A(1,2) = 1.0/h^2;
        for j = 2:m-1
            A(j,[j-1,j+1]) = [1.0/h^2, 1.0/h^2];
        end
        A(m,m-1) = 1.0/h^2;
    end

% Newton's method
printf('k        ||delta^k||_inf\n')
tt = h:h:T-h;         % for plotting only
UU = U0(:);           % force into column
for k = 0:K-1
    plot(tt,UU,'k'),  hold on
    m3 = floor(m/3);  % location to print k in figure
    text(tt(m3),1.2*UU(m3),sprintf('%d',k),'fontsize',16)
    delta = - Jacobian(UU) \ G(UU);    % solve  J(U) delta = - G(U)
    UU = UU + delta;
    printf('%d        %.4e\n',k,norm(delta,'inf'))
end
plot(tt,UU,'k')
xlabel t,  hold off,  axis tight
end  % function pendulum
