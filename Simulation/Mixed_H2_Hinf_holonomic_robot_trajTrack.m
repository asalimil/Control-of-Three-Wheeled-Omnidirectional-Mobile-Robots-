clc; clear; close all;

X = [3.5; 4.5]; Xdot = [0; 0]; Xddot = [0; 0];
T = 1000;
k = 1; lambda = 1;

t = 0; dt = 0.1;
ne = inf;

O_obs = [-1; 0]; R_obs = 0.5;

%
while ne > 0.05
    
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
    % w2 = [0; 0];
    w2 = [0.2*sin(t/2); 0.3*sin(t)];
    
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
    
    %     Y = sdpvar(4,4);
    %     Z = sdpvar(2,4,'full');
    %     W = sdpvar(6,6);
    %     beta_d = sdpvar(1);
    %     eta = 0.1;
    %     Const = [];
    %     Const = [Const; Y >= eta*eye(size(Y))];
    %     M1 = A*Y+Bu*Z+Y*A'+Z'*Bu'+ Bw*Bw';
    %     M2 = [Y (C1*Y+D12*Z)'
    %          (C1*Y+D12*Z) W];
    %     Const = [Const; M1 <= 0 ; M2 >=0 ; trace(W) <= beta_d];
    %     optimize(Const, beta_d);
    %     value(beta_d);
    %     gamma_d = sqrt(value(beta_d));
    %     H2_optimal_gain=value(gamma_d);
    %     Y = value(Y);
    %     eig_Y = eig(Y);
    %     Z = value(Z);
    %     F = Z*pinv(Y);
    
    beta1=sdpvar(1);
    beta2=sdpvar(1);
    X1=sdpvar(4);
    Y1=sdpvar(4);
    Z=sdpvar(6);
    An=sdpvar(4,4);
    Bn=sdpvar(4,8,'full');
    Cn=sdpvar(2,4,'full');
    Dn=sdpvar(2,8);
    Const=[];
    M1=[A*Y1+Y1*A'+Bu*Cn+Cn'*Bu'      (A'+An+(Bu*Dn*C2)')'               (Bw+Bu*Dn*D21);
        (A'+An+(Bu*Dn*C2)')           X1*A+A'*X1+Bn*C2+C2'*Bn'           (X1*Bw+Bn*D21);
        (Bw+Bu*Dn*D21)'                           (X1*Bw+Bn*D21)'                                     -eye(10)];
    
    M2=[Y1                                       eye(4)                (C1*Y1+D12*Cn)';
        eye(4)                                     X1                  (C1+D12*Dn*C2)';
        (C1*Y1+D12*Cn)       (C1+D12*Dn*C2)                                Z];
    
    Const=[Const;  (D11+D12*Dn*D21)==0; trace(Z) <= beta1];
    
    M3=[A*Y1+Y1*A'+Bu*Cn+Cn'*Bu'   (A'+An+(Bu*Dn*C2)')'   (Bw+Bu*Dn*D21)           (C1*Y1+D12*Cn)';
        (A'+An+(Bu*Dn*C2)')   X1*A+A'*X1+Bn*C2+C2'*Bn'     (X1*Bw+Bn*D21)         (C1+D12*Dn*C2)';
        (Bw+Bu*Dn*D21)'              (X1*Bw+Bn*D21)'               -beta2*eye(10)                  (D11+D12*Dn*D21)';
        (C1*Y1+D12*Cn)             (C1+D12*Dn*C2)              (D11+D12*Dn*D21)                             -eye(6)];
    
    Const=[Const; M1 <= 0];
    Const=[Const; M2 >= 0];
    Const=[Const; M3 <= 0];
    beta=beta1+beta2;
    optimize(Const,beta);
    beta=value(beta);
    value(beta1);
    value(beta2);
    gamma1=sqrt(value(beta1));
    gamma2=sqrt(value(beta2));
    H_2_optimal_gain=value(gamma1)
    H_infinity_optimal_gain=value(gamma2)
    
    %
    X1=value(X1);
    Y1=value(Y1);
    X2=eye(4)-X1*Y1;
    Y2=eye(4);
    An=value(An);
    Bn=value(Bn);
    Cn=value(Cn);
    Dn=value(Dn);
    % Calculating Ak2 Bk2 Ck2 Dk2 from An Bn Cn Dn X1 Y1 found in LMI
    M3 = pinv([X2 X1*Bu;zeros(2,4) eye(2)])* ... 
        ([An Bn; Cn Dn]-[X1*A*Y1 zeros(4,8); zeros(2,4) zeros(2,8)])* ... 
        pinv([Y2' zeros(4,8); C2*Y1 eye(8)]);
    Ak2=M3(1:4,1:4); Bk2=M3(1:4,5:12);
    Ck2=M3(5:6,1:4); Dk2=zeros(2,8);
    % Dk2=M3(7:9,7:9);
    % Calculating Ak Bk Ck Dk from Ak2 Bk2 Ck2 Dk2
    Dk=pinv(eye(2)+Dk2*D22)*Dk2;
    Ck=(eye(2)-Dk*D22)*Ck2;
    Bk=Bk2*(eye(8)-D22*Dk);
    Ak=Ak2-Bk*pinv(eye(8)-D22*Dk)*D22*Ck;
    F = [Ak Bk; Ck Dk];
    
    % ---
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
    drawnow;
    % axis([-5 5 -5 5]);
    axis equal
    %
    ne = norm(X(1:2));
    t = t + dt;
end
