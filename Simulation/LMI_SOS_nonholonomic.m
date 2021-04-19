clc; clear; close all;

% system
sdpvar ro beta
kr = 0.2; kq = 0.5;
f = [kr*ro*cos(beta)
     (kq-kr)*sin(beta)];
[V, Vc] = polynomial([ro beta],4);
F = Vc(1) == 0;
F = [F; sos(V - 0.00001*(ro^2 + beta^2))];
nablaV = jacobian(V,[ro beta]);
F = [F; sos(-nablaV*f)];
solvesos(F,[],[],[Vc])
