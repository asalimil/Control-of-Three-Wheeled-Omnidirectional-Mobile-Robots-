clc; clear; close all;

X = [3.5; 4.5]; Xdot = [0; 0]; Xddot = [0; 0];
T = 1000; k = 1; lambda = 1;
t = 0; dt = 0.1;
ne = inf;

while ne > 0.05
    
    O_obs = [0; -1]; R_obs = 0.25;
    
    
    A = [0         0   1   0
        0         0   0   1
        0         0   0   0
        0         0   0   0];
    m0 = 1;
    m = m0*(0.1 + exp(-1*t));
    B = [0 0
        0 0
        1/m0 0
        0 1/m0];
    C = eye(4);
    D = [0 0
        0 0
        0 0
        0 0];
    
    Bw=[zeros(4,4) B zeros(4,4)];
    Bu=B;
    C1=[-C;zeros(2,4)];
    C2=[zeros(4,4); C];
    D11=[eye(4) -D zeros(4,4); zeros(2,4) zeros(2,2) zeros(2,4)];
    D12=[-D;eye(2)];
    D21=[eye(4) zeros(4,2) zeros(4,4); zeros(4,4) D eye(4)];
    D22=[zeros(4,2);D];
    P=[A Bw Bu; C1 D11 D12; C2 D21 D22];
    
    % w2 = [0.5; 0.3];
    % w2 = [0.0; 0.0];
    % w1 = [0.2; 0.5];
    
    % w2 = [0.5*sin(t); 0.3*sin(0.5*t)];
    w2 = [0; 0];
    Xr = [cos(0.5*t); sin(0.5*t)];
    Xrdot = [-0.5*sin(0.5*t); 0.5*cos(0.5*t)];
    Xrddot = [-0.25*cos(0.5*t); -0.25*sin(0.5*t)];
    
    % LMI
    Y = sdpvar(4,4);
    W = sdpvar(2,4);
    gamma = sdpvar(1);
    eta = 0.1;
    Const = [];
    %
    % X = sdpvar(2,1);
    Xdot = sdpvar(2,1);
    h = (X - O_obs)'*(X - O_obs) - R_obs^2;
    hdot = 2*(X - O_obs)'*Xdot;
    Beta = -log(h/(1+h));
    alpha = 1;
    %
    Const = [Const; Y >= eta*eye(size(Y)); -hdot/(h + h^2) - (alpha/Beta) <= 0];
    M = [Y*A'+A*Y+W'*Bu'+Bu*W      Bw                 Y*C1'+W'*D12'
        Bw'                   -gamma*eye(10)          D11'
        C1*Y+D12*W                D11                -gamma*eye(6)];
    Const = [Const; M <= 0];
    optimize(Const, gamma);
    Y = value(Y); W = value(W);
    %
    % X = value(X);
    Xdot = value(Xdot);
    %
    F = W*pinv(Y)
    
    U_tilda = F*[X - Xr; Xdot - Xrdot];
    
    % S = (Xdot - Xrdot) + lambda*(X - Xr);
    % U = -k*S + Xrddot;
    U = [0; 0] + Xrddot;
    Xddot = U + U_tilda + w2;
    Xdot = Xdot + Xddot*dt;
    X = X + Xdot*dt;
    %
    plot(X(1),X(2),'*r'); grid on; hold on;
    plot(Xr(1),Xr(2),'*g'); grid on; hold on;
    boundarycolor = [102, 102, 0]; fill_obs = 1;
    plot_circle(O_obs,R_obs,boundarycolor,fill_obs);

    drawnow;
    % axis([-5 5 -5 5]);
    axis equal
    %
    ne = norm(X(1:2));
    t = t + dt;
end
