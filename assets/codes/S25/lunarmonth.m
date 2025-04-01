% LUNARMONTH  Solve Earth-Moon problem using N=960 RK4, and removing
% center-of-mass, to compute the length of a lunar month in days.

m1 = 5.972e24;  m2 = 7.348e22;
eta = [0; 0; 3.944e8; 0;          % positions x1 y1 x2 y2
       0; 0;       0; 1.022e3];   % velocities v1 w1 v2 w2
t0 = 0;  secperday = 24 * 60 * 60;  tf = 40.0 * secperday;

% for center-of-mass straight-line motion
cx0 = (eta(1)*m1 + eta(3)*m2)/(m1+m2);
cy0 = (eta(2)*m1 + eta(4)*m2)/(m1+m2);
cv  = (eta(5)*m1 + eta(7)*m2)/(m1+m2);
cw  = (eta(6)*m1 + eta(8)*m2)/(m1+m2);

% hourly RK4
[tt,uu] = rk4(@fearthmoon,eta,t0,tf,960);
cxy = [cx0 + cv * tt; cy0 + cw * tt];   % center-of-mass pos.
xy = uu(1:4,:) - [cxy; cxy];  % megameters

% distance of moon from its initial position
d = sqrt((xy(3,:)-xy(3,1)).^2 + (xy(4,:)-xy(4,1)).^2);

% find minimizing index j
[z,j] = min(d(5:end));   % don't look at first 4 hours
(j+4) / 24               % in days
