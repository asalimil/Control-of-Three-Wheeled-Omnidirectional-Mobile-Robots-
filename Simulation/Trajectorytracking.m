clc; clear; close all;

xmin=0.5; % setting the figure limits
xmax=4;
ymin=-1;
ymax=5;


X = [1; 1];
Xdot = [0; 0];
U = [0; 0];
t = 0; dt = 0.01; T = 200;
k = 1; lambda = 5; r = 1;
a1 = 2; a2 = 0.5; a3 = 7;
cplot = 1; ne = inf;

% ne_max = ne;
while (ne > 0.05) || (t < T)
    
    % disturbance
    d = [0.1*sin(t/2) ; 0.3*sin(t/4)];
    n = [0.3*sin(2*t) ; 0.2*sin(3*t); 0.5*sin(4*t) ; 0.25*sin(2*t)];
    
    %     d = [0; 0; 0; 0];
    %     n = [0; 0; 0; 0];
    
    % Desired trajectory
    Xr = [2+cos(t); 2+sin(t)];
    Xrdot = [-sin(t); cos(t)];
    Xrddot = [-cos(t); -sin(t)];
    
    % Desired trajectory
    %     Xr = [a1*sin(t/a3); a2*sin(t/(2*a3))];
    %     Xrdot = [(a1/a3)*cos(t/a3); (a2/(2*a3))*cos(t/(2*a3))];
    %     Xrddot = [-(a1/(a3^2))*sin(t/a3); -(a2/(4*(a3^2)))*sin(t/(2*a3))];
    
    e = (X+n(1:2)) - Xr;
    edot = (Xdot+n(3:4)) - Xrdot;
    S = edot + lambda*e;
    U = -k*S + Xrddot;
    %
    
    Xddot = U + d;
    Xdot = Xdot + Xddot*dt;
    X = X + Xdot*dt;
    %
    %     ne = norm(e)
    %     drawnow;
    
    % Plot robot
    subplot(1,1,1)
    if (cplot == 1)
        c = 1;
        draw_robot_omni(X(1),X(2),0,c);
    elseif (mod(cplot,50) == 0) && (ne > 0.01)
        c = 0;
        draw_robot_omni(X(1),X(2),0,c);
    end
    if ne > 0.05
        plot(X(1),X(2),'.b'); hold on; grid on;
    else
        plot(X(1),X(2),'.r'); hold on; grid on;
    end
    
    % Plot X & Xr
    plot(X(1),X(2),'.r','LineWidth',1); grid on; drawnow; hold on;
    plot(Xr(1),Xr(2),'.g','LineWidth',1); grid on; drawnow; hold on;
    xlabel('x (m)'); ylabel('y (m)');
    axis([xmin xmax ymin ymax]) % setting the figure limits
    axis square
    % axis equal
    
    ne = norm(e);
    nedot = norm(edot);
    
    
    %     % Plot errors
    %     subplot(1,2,2)
    %     plot(t,ne,'*b','LineWidth',1); hold on; grid on; drawnow;
    %     xlabel('t (s)'); ylabel('norm of error (m)');
    %     axis square
    
    % plot(t,ne,'*b'); grid on; drawnow; hold on;
    %     plot(t,nedot,'*r'); grid on; drawnow; hold on;
    cplot = cplot + 1;
    t = t + dt;
    
end