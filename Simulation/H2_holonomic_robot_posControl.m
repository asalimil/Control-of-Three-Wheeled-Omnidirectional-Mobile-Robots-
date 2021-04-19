clc; clear; close all;

X = [3.5; 2]; Xdot = [0; 0]; Xddot = [0; 0];
T = 1000;
k = 1; lambda = 1;

t = 0; dt = 0.1;
ne = inf;
O_obs = [-1; 0]; R_obs = 0.5;

%
while ne > 0.05
    
    d = [0.3*sin(t) ; 0.2*sin(t/2)];
    n = [0.1; 0.3];
    
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
    
    % Regulator
    Bw=[zeros(4) B zeros(4)];
    Bu=B;
    C1=[-C; zeros(2,4)];
    C2=C;
    D11=[eye() D zeros(4);zeros(2,2) zeros(2,4) zeros(2,4)];
    D12=[-D;eye(2)];
    D21=[eye() D eye(4)];
    D22=D;
    P=[A Bw Bu
       C1 D11 D12
       C2 D21 D22];

   
   % w2 = [0.5; 0.3];
    % w2 = [0.0; 0.0];
    % w1 = [0.2; 0.5];
    
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
    
    % reference trajectory
    % Xr = [2+cos(t); 2+sin(t)];
    % Xrdot = [-sin(t); cos(t)];
    % Xrddot = [-cos(t); -sin(t)];
    Xr = [0; 0];
    Xrdot = [0; 0];
    Xrddot = [0; 0];
    
    %     % LMI
    %     Y = sdpvar(4,4);
    %     W = sdpvar(2,4);
    %     gamma = sdpvar(1);
    %     eta = 0.1;
    %     Const = [];
    %     %
    %     Const = [Const; Y >= eta*eye(size(Y))];
    %     M = [Y*A'+A*Y+W'*Bu'+Bu*W      Bw                 Y*C1'+W'*D12'
    %         Bw'                   -gamma*eye(10)          D11'
    %         C1*Y+D12*W                D11                -gamma*eye(6)];
    %     Const = [Const; M <= 0];
    %     optimize(Const, gamma);
    %     Y = value(Y); W = value(W);
    %     %
    %     F = W*pinv(Y)
    
    Y = sdpvar(4,4);
    Z = sdpvar(2,4,'full');
    W = sdpvar(6,6);
    beta_d = sdpvar(1);
    eta = 0.1;
    Const = [];
    Const = [Const; Y >= eta*eye(size(Y))];
    M1 = A*Y+Bu*Z+Y*A'+Z'*Bu'+ Bw*Bw';
    M2 = [Y (C1*Y+D12*Z)'
        (C1*Y+D12*Z) W];
    Const = [Const; M1 <= 0 ; M2 >=0 ; trace(W) <= beta_d];
    optimize(Const, beta_d);
    value(beta_d);
    gamma_d = sqrt(value(beta_d));
    H2_optimal_gain=value(gamma_d);
    Y = value(Y);
    eig_Y = eig(Y);
    Z = value(Z);
    F = Z*pinv(Y);
    
    
    U_tilda = F*[(X + n) - Xr; Xdot - Xrdot];
    
    % S = (Xdot - Xrdot) + lambda*(X - Xr);
    % U = -k*S + Xrddot;
    U = [0; 0] + Xrddot + d;
    Xddot = U + U_tilda;
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
