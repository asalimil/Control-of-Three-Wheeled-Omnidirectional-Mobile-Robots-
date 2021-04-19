% LMI for Quadcopter control

clc; clear; close all;

X = [3.5; 4.5]; Xdot = [0; 0]; Xddot = [0; 0];
T = 1000;
k = 1; lambda = 1;

t = 0; dt = 0.1;
ne = inf;

O_obs = [-1; 0]; R_obs = 0.5;

%
while ne > 0.05
    
%     A = [0         0   1   0
%         0         0   0   1
%         0         0   0   0
%         0         0   0   0];
%     m0 = 1;
%     m = m0*(0.1 + exp(-1*t));
%     B = [0 0
%         0 0
%         1/m0 0
%         0 1/m0];
%     C = eye(4);
%     D = [0 0
%         0 0
%         0 0
%         0 0];
    g = 9.8;
    A = [0 0 0 1 0 0 0 0 0 0 0 0
         0 0 0 0 1 0 0 0 0 0 0 0 
         0 0 0 0 0 1 0 0 0 0 0 0 
         0 0 0 0 0 0 0 0 0 0 0 0 
         0 0 0 0 0 0 0 0 0 0 0 0 
         0 0 0 0 0 0 0 0 0 0 0 0 
         0 -g 0 0 0 0 0 0 0 0 0 0 
         g 0 0 0 0 0 0 0 0 0 0 0 
         0 0 0 0 0 0 0 0 0 0 0 0 
         0 0 0 0 0 0 1 0 0 0 0 0 
         0 0 0 0 0 0 0 1 0 0 0 0 
         0 0 0 0 0 0 0 0 1 0 0 0 
        ];
   Ix = 1; Iy = 1; Iz = 1; m = 1;
   B = [0 0 0 0
        0 0 0 0
        0 0 0 0
        0 1/Ix 0 0 
        0 0  1/Iy 0
        0 0   0  1/Iz
        0 0 0 0
        0 0 0 0
        1/m 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0];
    % fwx fwy fwz tau_wx tau_wy tau_wz
    d = [0 0 0 0 0 0]';
    w1 = [0 0 0 0 0 0
         0 0 0 0 0 0
         0 0 0 0 0 0
         0 0 0 1/Ix 0 0 
         0 0 0 0 1/Iy 0
         0 0 0 0 0 1/Iz
         1/m 0 0 0 0 0
         0 1/m 0 0 0 0
         0 0 1/m 0 0 0
         0 0 0 0 0 0
         0 0 0 0 0 0
         0 0 0 0 0 0];
    C = eye(12);
    D = zeros(12,4);
    % System
    % xdot = Ax + Bu +w1d
    
    
    Bw=[zeros(12,3) B zeros(12,3)];
    Bu=B;
    C1=[-C;zeros(4,12)];
    C2=[zeros(3,12); C];
    D11=[eye(3) -D zeros(3,3); zeros(4,3) zeros(4,4) zeros(4,3)];
    D12=[-D;eye(4)];
    D21=[eye(3) zeros(3,4) zeros(3,3); zeros(3,3) D eye(3)];
    D22=[zeros(3,4);D];
    P=[A Bw Bu 
       C1 D11 D12
       C2 D21 D22];
       
    % w2 = [0.5*sin(t); 0.3*sin(0.5*t)];
    % w2 = [0; 0];
    % w2 = [0.2*sin(t/2); 0.3*sin(t)];
    
    %
    %     Xr = [cos(0.5*t); sin(0.5*t)];
    %     Xrdot = [-0.5*sin(0.5*t); 0.5*cos(0.5*t)];
    %     Xrddot = [-0.25*cos(0.5*t); -0.25*sin(0.5*t)];
    %     Xr = [2*sin(t/7); 0.5*sin(t/(14))];
    %     Xrdot = [(2/7)*cos(t/7); (0.5/14)*cos(t/14)];
    %     Xrddot = [-(2/49)*sin(t/7); -(0.5/(14*14))*sin(t/14)];
    Xr = [2+cos(t); 2+sin(t)];
    Xrdot = [-sin(t); cos(t)];
    Xrddot = [-cos(t); -sin(t)];
    
    % LMI
    Y = sdpvar(12,12);
    W = sdpvar(4,12);
    gamma = sdpvar(1);
    eta = 0.1;
    Const = [];
    %
    Const = [Const; Y >= eta*eye(size(Y))];
    M = [Y*A'+A*Y+W'*Bu'+Bu*W      Bw                 Y*C1'+W'*D12'
        Bw'                   -gamma*eye(10)          D11'
        C1*Y+D12*W                D11                -gamma*eye(7)];
    Const = [Const; M <= 0];
    optimize(Const, gamma);
    Y = value(Y); W = value(W);
    %
    F = W*pinv(Y)
    
    U_tilda = F*[X - Xr; Xdot - Xrdot];
    
    % S = (Xdot - Xrdot) + lambda*(X - Xr);
    % U = -k*S + Xrddot;
    U = [0; 0] + Xrddot + U_tilda;
    % Xddot = U + w2;
    Xdot = A*X + B*U + w1*d;
    Xdot = Xdot + Xddot*dt;
    X = X + Xdot*dt;
    %
    plot(X(1),X(2),'*r'); grid on; hold on;
    plot(Xr(1),Xr(2),'*g'); grid on; hold on;
    drawnow;
    % axis([-5 5 -5 5]);
    axis equal
    %
    ne = norm(X(1:2));
    t = t + dt;
end
