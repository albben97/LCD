%% Q Ax1
clc; clear all;
syms i s phi1 phi2 phi3 omega1 omega2 omega3 real
syms R Ke Kt J1 J2 b d1 d2 real
x = [phi1 phi2 phi3 omega1 omega2];
R = 1005.03; 
%syms R real positive
Ke = 0.1; Kt = 0.1; 
J1 = 1*10^-5; J2 = 4*10^-5;
b = 2 * 10^-3; 
d1 = 20; d2 = 2;

A = [0 0 0 1 0;
   0 0 0 0 1;    
   0 d2/b -d2/b 0 0; 
   -d1/J1 d1/J1 0 -(Ke*Kt)/(R*J1) 0; 
   d1/J2 -(d2+d1)/J2 d2/J2 0 0];

B = [0 0; 0 0; 0 1/b; Kt/(R*J1) 0; 0 0];
C1 = [0 1 0 0 0; 0 0 0 0 1];
D1 = [0 0; 0 0];
C2 = [0 0 0 -Ke/R 0; 0 d2/b -d2/b 0 0];
D2 = [1/R 0; 0 1/b];
C = C2;
D = D2;

%observability 
%R > 0.0054 gives rank 4 for C1
%R > 0.0078 gives rank 4 for C2 
obsv1 = [C; C*A; C*A*A; C*A*A*A; C*A*A*A*A];
obsv2 = rref(obsv1)
obsvctb = obsv(A, C);
rref(obsvctb)
rank(obsvctb);
%obsv = rref(obsv)
%obsv = simplify(obsv)

%controllability
ctrb = [B A*B A*A*B A*A*A*B A*A*A*A*B];
ctrb1 = rref(ctrb);
%ctrb = simplify(cntr)

%sys = ss(A,B,C2,D2);
%obsv_sys = obsv(sys)
rank_obsv_sys = rank(obsv2);

%ctrb_sys = ctrb(sys);
%rank_ctrb_sys = rank(ctrb1)
%% Discrete (d) & e))
clc
% Sampling interval (s)
Ts = 10^(-3);

% Discrete A-matrix
Ad = expm(A*Ts)

% Discrete B-matrix, constructed using ZOH
syms t real positive
fun = expm(A.*t);
%Bd = int(fun,[0 Ts])*B;
lower = expm(A*0)
upper = expm(A*Ts)

inv(A)*(upper - lower)*B

Bd = inv(A)*(Ad-eye(5))*B

%syms s t
%matrix = inv(eye(size(A))*s-A);
%exp_At = vpa(ilaplace(matrix));

%Bd = int(exp_At,t,0,1e-3)*B
%% Discrete second attempt
clc
% Sampling interval (s)
Ts = 10^(-3);

% Discrete A-matrix
Ad = expm(A*Ts)

% Discrete B-matrix, constructed using ZOH
syms t real positive
fun = expm(A.*t);
%Bd = int(fun,[0 Ts])*B;
%lower = expm(A*0)
%upper = expm(A*Ts)
f = @(t)expm(A*t)*B;
Bdalbert = integral(f,0,Ts,'ArrayValued',true)
inv(A)*(upper - lower)*B

Bd = inv(A)*(Ad-eye(5))*B

%syms s t
%matrix = inv(eye(size(A))*s-A);
%exp_At = vpa(ilaplace(matrix));
%%
Ad = expm(A*Ts)
Astar = adjoint(A)
C= C2
obsvd = [C; C*Ad; C*Ad*Ad; C*Ad*Ad*Ad; C*Ad*Ad*Ad*Ad];
obsvd = rref(obsvd)