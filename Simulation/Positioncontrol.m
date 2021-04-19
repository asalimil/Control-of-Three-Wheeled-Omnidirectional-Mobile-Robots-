clc; clear; close all;

% Regulator problem
xmin=-0.5; % setting the figure limits
xmax=2;
ymin=-0.5;
ymax=2;

X = [1.5; 1];
Xdot = [0; 0];
U = [0; 0];
t = 0; dt = 0.05; T = 10;
lambda = 0.5; k = 2;
r = 0.5;
ne = inf;

rgb = [46,139,87];
plot(0,0,'*','Color',rgb/255,'LineWidth',2); hold on

plotted = 0;
cplot = 1;
while (ne > 0.01) || (t < T)
    
    d = [0.3*sin(t) ; 0.2*sin(t/2)];
    n = [0.1; 0.3];
    
    
    Xr = [0; 0];
    Xrdot = [0; 0];
    Xrddot = [0; 0];
    %
    e = (X + n ) - Xr;
    edot = Xdot - Xrdot;
    S = edot + lambda*e;
    U = -k*S + Xrddot;
    %
    Xddot = U + d;
    Xdot = Xdot + Xddot*dt;
    X = X + Xdot*dt;
    %
    subplot(1,1,1)
    if (cplot == 1)
        plotted = plotted + 1;
        c = 1;
        draw_robot_omni(X(1),X(2),0,c);
    elseif (mod(cplot,10) == 0) && (ne > 0.01) && plotted < 12
        plotted = plotted + 1;
        c = 0;
        draw_robot_omni(X(1),X(2),0,c);
    end
    plot(X(1),X(2),'.r','LineWidth',3); grid on; drawnow;
    rgb = [46,139,87];
    plot(0,0,'*','Color',rgb/255,'LineWidth',2); hold on
    xlabel('x (m)'); ylabel('y (m)');
    axis([xmin xmax ymin ymax]) % setting the figure limits
    axis square
    % axis equal
    %
    ne = norm(X);
%     subplot(2,2,3)
%     plot(t,ne,'.b','LineWidth',3); hold on; grid on; drawnow;
%     xlabel('t (s)'); ylabel('norm of error (m)');
%     axis square
    cplot = cplot + 1;
    t = t + dt;
    
end


