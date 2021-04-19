clc; clear; close all;

X = [3.5; 4.5]; 
Xdot = [0; 0];
U = [0; 0];
t = 0; dt = 0.1; T = 500;
k = 1; lambda = 0.5; r = 1;
a1 = 2; a2 = 0.5; a3 = 7;

while t < T
    
    % disturbance
    % w1 = [0.2*sin(t/2); 0.3*sin(t)];
    w1 = [0.2; 0.3];
    
    % w.o disturbance 
    % w1 = [0; 0];
    
    t = t + dt;
    % Desired trajectory
    Xr = [2+cos(t); 2+sin(t)];
    Xrdot = [-sin(t); cos(t)];
    Xrddot = [-cos(t); -sin(t)];
    
    % Desired trajectory
%     Xr = [a1*sin(t/a3); a2*sin(t/(2*a3))];
%     Xrdot = [(a1/a3)*cos(t/a3); (a2/(2*a3))*cos(t/(2*a3))];
%     Xrddot = [-(a1/(a3^2))*sin(t/a3); -(a2/(4*(a3^2)))*sin(t/(2*a3))];
    
    e = X - Xr;
    edot = Xdot - Xrdot;
    S = edot + lambda*e;
    U = -k*S + Xrddot;
    %
    
    Xddot = U + w1;
    Xdot = Xdot + Xddot*dt;
    X = X + Xdot*dt;
    %
    ne = norm(e)
    % if ne > 0.05
        plot(X(1),X(2),'*r'); hold on; grid on; axis equal;
    % else
    %     plot(X(1),X(2),'*r'); hold on; grid on; axis equal;
    % end
    plot(Xr(1),Xr(2),'*g'); hold on; grid on; axis equal;
    drawnow;
end