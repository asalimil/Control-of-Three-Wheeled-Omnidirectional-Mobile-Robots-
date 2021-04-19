clc; clear; close all;

Z = [3.5; 4.5]; Zdot = [0; 0]; Zddot = [0; 0];
T = 1000;
k = 1; lambda = 1;

A = [0         0   1   0
    0         0   0   1
    -k*lambda 0   -k  0
    0         -k*lambda 0 -k];
B = [0 0
    0 0
    1 0
    0 1];
C = eye(4);
D = [0 0
    0 0
    0 0
    0 0];

Bw = [B zeros(4,4)];
Bu = B;
C1 = [C; zeros(2,4)];
C2 = C;
D11 = [D zeros(4,4); zeros(2,2) zeros(2,4)];
D12 = [D; eye(2)];
D21 = [D eye(4)];
D22 = D;
P = [A  Bw  Bu
    C1 D11 D12
    C2 D21 D22];
t = 0; dt = 0.1;
ne = inf;
%
while ne > 0.05
    Zr = [0; 0]; 
    Zrdot = [0; 0]; 
    Zrddot = [0; 0];
    % LMI
    Y = sdpvar(4,4);
    W = sdpvar(2,4);
    gamma = sdpvar(1);
    eta = 0.000001;
    %
    Const = [];
    Const = [Const; Y >= eta*eye(size(Y))];
    M = [Y*A'+A*Y+W'*Bu'+Bu*W      Bw                 Y*C1'+W'*D12'
       Bw'                   -gamma*eye(6)          D11'
       C1*Y+D12*W                D11                -gamma*eye(6)];
    Const = [Const; M <= 0];
    optimize(Const, gamma);
    Y = value(Y); W = value(W);
    %
    F = W*pinv(Y);
    U_tilda = F*[Z; Zdot];
    %
    S = Zdot + lambda*Z;
    U = -k*S + Zrddot;
    Zddot = U + U_tilda;
    Zdot = Zdot + Zddot*dt;
    Z = Z + Zdot*dt;
    %
    plot(Z(1),Z(2),'.r'); grid on; hold on; drawnow;
    axis([-5 5 -5 5]);
    %
    ne = norm(Z(1:2));
end
