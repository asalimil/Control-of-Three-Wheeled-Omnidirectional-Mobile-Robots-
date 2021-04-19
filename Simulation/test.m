clc; clear; close all;

X = [3.5; 4.5]; Xdot = [0; 0]; Xddot = [0; 0];
k = 1; lambda = 1;

t = 0; dt = 0.1;
ne = inf;

while ne > 0.05
    t = t + dt;
    Xr = [cos(t); sin(t)]; 
    Xrdot = [-sin(t); cos(t)]; 
    Xrddot = [-cos(t); -sin(t)];
    S = (Xdot - Xrdot) + lambda*(X - Xr);
    U = -k*S + Xrddot;
    Xddot = U;
    Xdot = Xdot + Xddot*dt;
    X = X + Xdot*dt;
    plot(X(1),X(2),'*r'); grid on; hold on; drawnow;
    ne = norm(X(1:2));
    axis equal
end
