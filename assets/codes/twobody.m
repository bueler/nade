% TWOBODY  Solve Earth-Moon problem by comparing forward Euler, RK4,
% and ODE45 solutions of a system of 8 first-order ODEs.  Removes
% the center-of-mass motion before plotting.
% Requires: FEARTHMOON

eta = [0; 0; 3.944e8; 0;          % positions x1 y1 x2 y2
       0; 0;       0; 1.022e3];   % velocities v1 w1 v2 w2
t0 = 0;
secperday = 24 * 60 * 60;
tf = 40.0 * secperday;

% for center-of-mass straight-line motion
% (see https://en.wikipedia.org/wiki/Two-body_problem)
m1 = 5.972e24;                    % mass of Earth in kg
m2 = 7.348e22;                    % mass of Moon in kg
cx0 = (eta(1)*m1 + eta(3)*m2)/(m1+m2);
cy0 = (eta(2)*m1 + eta(4)*m2)/(m1+m2);
cv  = (eta(5)*m1 + eta(7)*m2)/(m1+m2);
cw  = (eta(6)*m1 + eta(8)*m2)/(m1+m2);

N = [40 960];
styles = {{'g.','r.','k*'};
          {'b.-','m.-','k*'}};

% k = 1 is daily figure, k = 2 is hourly figure
for k = 1:2
    figure(k),  hold on
    for m = 1:3  % m=1,2,3 is Euler,RK4,ode45
        if m==1
            [tt,uu] = feuler(@fearthmoon,eta,t0,tf,N(k));
        elseif m==2
            [tt,uu] = rk4(@fearthmoon,eta,t0,tf,N(k));
        elseif m==3
            [tt,uu] = ode45(@fearthmoon,[t0,tf],eta);
            % optionally with higher accuracy:
            %[tt,uu] = ode45(@fearthmoon,[t0,tf],eta,...
            %                odeset('RelTol',1.0e-10,'AbsTol',1.0e-10));
            tt = tt';  uu = uu';
        end
        cxy = [cx0 + cv * tt; cy0 + cw * tt];   % center-of-mass pos.
        xy = (uu(1:4,:) - [cxy; cxy]) / 1.0e6;  % megameters
        if m < 3
            plot(xy(1,:),xy(2,:),styles{1}{m},...
                 xy(3,:),xy(4,:),styles{2}{m})
        else
            plot(xy(3,:),xy(4,:),styles{2}{m})
        end
    end
    hold off
    xlabel('x  (Mm)'),  ylabel('y  (Mm)')
    legend('FE (earth)','FE (moon)','RK4 (earth)','RK4 (moon)',...
           'ode45 (moon)','Location','SouthWest');
    axis([-900 700 -900 700]),  grid on,  axis square
end
figure(1),  print -dpng dailytwobody.png
figure(2),  print -dpng hourlytwobody.png
