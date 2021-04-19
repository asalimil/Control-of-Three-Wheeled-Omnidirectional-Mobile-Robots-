clc; clear; close all;

X = [3.5; 4.5]; Xdot = [0; 0]; Xddot = [0; 0];
T = 1000;
k = 1; lambda = 1;


t = 0; dt = 0.05;
ne = inf;

O_obs1 = [0; -1]; R_obs1 = 0.25;
% O_obs2 = [1; 0]; R_obs2 = 0.25;

%
% while ne > 0.05
while t < T
    
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
    Xdot = sdpvar(2,1);
    H1 = 2*(X - O_obs1)'*Xdot;
    % H2 = 2*(X - O_obs2)'*Xdot;
    % H1 = (X - O_obs1)'*(X - O_obs1) - R_obs1^2;
    % H2 = (X - O_obs2)'*(X - O_obs2) -  R_obs2^2;
    %
    % Const = [Const; Y > eta*eye(size(Y)); H1 > zeros(size(H1)); H2 > zeros(size(H2))];
    Const = [Const; Y > eta*eye(size(Y)); H1 > zeros(size(H1))];
    M = [Y*A'+A*Y+W'*Bu'+Bu*W      Bw                 Y*C1'+W'*D12'
        Bw'                   -gamma*eye(10)          D11'
        C1*Y+D12*W                D11                -gamma*eye(6)];
    Const = [Const; M <= 0];
    optimize(Const, gamma);
    Y = value(Y); W = value(W);
    %
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
    plot(X(1),X(2),'.r'); grid on; hold on;
    plot(Xr(1),Xr(2),'.g'); grid on; hold on;
    boundarycolor = [102, 102, 0]; fill_obs = 1;
    plot_circle(O_obs1,R_obs1,boundarycolor,fill_obs);
%     boundarycolor = [102, 102, 0]; fill_obs = 1;
%     plot_circle(O_obs2,R_obs2,boundarycolor,fill_obs);

    drawnow;
    % axis([-5 5 -5 5]);
    axis equal
    %
    ne = norm(X(1:2));
    t = t + dt;
end
