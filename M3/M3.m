% Assignment M1
beep off; clc; clear all;

% Variables
syms phi_1 phi_2 phi_3 w_1 w_2 real

% Parameter Values
R = 1;  %[Ohm]
K_E = 10^(-1);  %[V*s/rad] 
K_T = 10^(-1);  %[Nm/A]
J_1 = 10^(-5);  %[kg*m^2]
J_2 = 4*10^(-5);    %[kg*m^2]
b = 2*10^(-3);  %[Nm*s]
D_1 = 20;   %[Nm/rad]
D_2 = 2;    %[Nm/rad]

% A-matrix from 1. b)
A = [0 0 0 1 0;
     0 0 0 0 1;
     0 D_2/b -D_2/b 0 0;
     -D_1/J_1 D_1/J_1 0 -K_E*K_T/(R*J_1) 0;
     D_1/J_2 -(D_2+D_1)/J_2 D_2/J_2 0 0];

% B-matrix from 1. b)
B = [0 0 0 K_T/(R*J_1) 0;
     0 0 1/b 0 0]';

% C-matrixes from 1. c)
C_1 = [0 1 0 0 0;
       0 0 0 0 1];
   
C_2 = [0 0 0 -K_E/R 0;
       0 D_2/b -D_2/b 0 0];
   
% D-matrixes from 1. c)
D_1 = [0 0;
       0 0];

D_2 = [1/R 0;
       0 1/b];
   
% State space vectors
x = [phi_1 phi_2 phi_3 w_1 w_2]';

% Eigenvalues for A
eig_A = eig(A)

% State-space of system
sys_1 = ss(A,B,C_1,D_1);
sys_2 = ss(A,B,C_2,D_2);

% Transfer functions
s = tf('s');
G_1 = C_1*inv(s*eye(5)-A)*B+D_1;
G_2 = C_2*inv(s*eye(5)-A)*B+D_2;

G_2_zpk = zpk(G_2)
poles = pole(G_2)
tzeros = tzero(G_2)

% Inputs for simulink
C = C_1;
D = D_1;
K = 1;

%% Assignment M2
beep off; clc; clear all;

% Variables
syms phi_1 phi_2 phi_3 w_1 w_2 real

syms R K_E K_T J_1 J_2 b D_1 D_2 real
%syms D_1 real
% Parameter Values
R = 1;  %[Ohm]
K_E = 10^(-1);  %[V*s/rad] 
K_T = 10^(-1);  %[Nm/A]
J_1 = 10^(-5);  %[kg*m^2]
J_2 = 4*10^(-5);    %[kg*m^2]
b = 2*10^(-3);  %[Nm*s]
D_1 = 20;   %[Nm/rad]
D_2 = 2;    %[Nm/rad]

% A-matrix from 1. b)
A = [0 0 0 1 0;
     0 0 0 0 1;
     0 D_2/b -D_2/b 0 0;
     -D_1/J_1 D_1/J_1 0 -K_E*K_T/(R*J_1) 0;
     D_1/J_2 -(D_2+D_1)/J_2 D_2/J_2 0 0];

% B-matrix from 1. b)
B = [0 0 0 K_T/(R*J_1) 0;
     0 0 1/b 0 0]';

% C-matrixes from 1. c)
C_1 = [0 1 0 0 0;
       0 0 0 0 1];
   
C_2 = [0 0 0 -K_E/R 0;
       0 D_2/b -D_2/b 0 0];
   
% D-matrixes from 1. c)
D1 = [0 0;
       0 0];

D2 = [1/R 0;
       0 1/b];

% State space vectors
x = [phi_1 phi_2 phi_3 w_1 w_2]';

% Eigenvalues for A
eig_A = eig(A);

% Observability matrices, both reduced matrices are displayed as well as
% their ranks. Rank = 4 means full rank in this case
%obsv1 = obsv(A,C_1);
%obsv1_reduced = rref(obsv1);
%obsv1_rank = rank(obsv1);

%obsv2 = obsv(A,C_2);
%obsv2_reduced = rref(obsv2);
%obsv2_rank = rank(obsv2);

% Controllability matrices, rank = 5 means full rank.
%ctrb = ctrb(A,B);
%ctrb_reduced = rref(ctrb);
%ctrb_rank = rank(ctrb);

% Transfer functions
%s = tf('s');
%G_1 = C_1*inv(s*eye(5)-A)*B+D1;
%G_2 = C_2*inv(s*eye(5)-A)*B+D2;

% PBH-test
%ctrb_PBH = [A-eig(A).*eye(size(A)) B];
%rank_ctrb_PBH = rank(ctrb_PBH);

%obsv1_PBH = [A-eig(A).*eye(size(A)); C_1];
%rank_obsv1_PBH = rank(obsv1_PBH);

%obsv2_PBH = [A-eig(A).*eye(size(A)); C_2];
%rank_obsv2_PBH = rank(obsv2_PBH);

%% Discrete (d) & e))
%clc
% Sampling interval (s)
Ts = 10^(-3);   

% Discrete A-matrix
Ad = expm(A*Ts)

% Discrete B-matrix, constructed using ZOH
Bd = inv(A)*(Ad-eye(5))*B

% Using integral, as suggested by Albert
f = @(t) expm(A*t)*B
Bd_int = integral(f,0,Ts,'ArrayValued',true)
Bd = Bd_int
%% f)
% Check stability of Ad using eigenvalues. If all eigen values are within
% unit circle the system is stable.

Ad_eig = eig(Ad);

% Observability matrices, both reduced matrices are displayed as well as
% their ranks. Rank = 4 means full rank in this case
%obsv1_d = obsv(Ad,C_1)
%obsv1_reduced_d = rref(obsv1_d)
%obsv1_rank_d = rank(obsv1_d)

%obsv2_d = obsv(Ad,C_2);
%obsv2_reduced_d = rref(obsv2_d)
%obsv2_rank_d = rank(obsv2_d)

% Controllability matrices, rank = 5 means full rank.
%ctrb_d = ctrb(Ad,Bd)
%ctrb_reduced_d = rref(ctrb_d)
%ctrb_rank_d = rank(ctrb_d)

%% Assingment M3
clc
% a)
PoR = 99.7;     % Percentage of Realization
mu = 0;
C = C_1;
syms sigma x real

norm_pdf = 1/(sqrt(2*pi*sigma))*exp(-1/2*(x-mu)^2/sigma);

% Voltage noises
Phi_v = int(norm_pdf, x, -0.3, 0.3);
sigma_val_v = double(solve(Phi_v==PoR/100,sigma));

% Torque noises
Phi_T = int(norm_pdf, x, -0.1, 0.1);
sigma_val_T = double(solve(Phi_T==PoR/100,sigma));

% measurement noises
Phi_T = int(norm_pdf, x, -0.02, 0.02);
sigma_val_noise1 = double(solve(Phi_T==PoR/100,sigma));

Phi_T = int(norm_pdf, x, -0.01, 0.01);
sigma_val_noise2 = double(solve(Phi_T==PoR/100,sigma));
R1 = sigma_val_v;
R2 = sigma_val_T;
QN = [R1 0; 0 R2];

RN = [sigma_val_noise1 0; 0 sigma_val_noise2];
NN = [0 0; 0 0];
%plant=ss(sysc.a,[sysc.b G],sysc.c,[sysc.d zeros(2)]);

h = 0.001
plant = ss(Ad, [Bd Bd], C_1, [D1 zeros(2)], h);

[filter,L_CT,P] = kalman(plant,QN,RN,NN);
%[K,P,Z,E] = dlqe(Ad,Bd,C,QN,RN)
eig(filter.a);
check =eig(Ad-L_CT*C)
N = Bd;
%
%A = [A(1,:) 0; A(2,:) 0; A(3,:) 0; A(4,:) 0; A(5,:) 0; C_1(1,:) 0];
%B = [B(1,:); B(2,:); B(3,:); B(4,:); B(5,:); 0 0];
%C = [C_1(1,:) 0; C_1(2,:) 0];
%C = [0 0 0 0 1 0];
%
C = [0 1 0 0 0; 0 0 0 0 1];
C_1 = C;
D = [0; 0; 0];
%Qx = diag([1/(1000^2) 1/(1000^2) 1/(1000^2) 10 10 0.001])
Qx = diag([1 10 1 1  1/1000 1/1000]) %phi1 phi2 phi3 w1 w2 int
Qu = diag([2 2]) * 1000000;
%K_CT = lqr(Ad,Bd,Qx,Qu)

%[K_CT,S,e] = dlqr(Ad,Bd,Qx,Qu)
%L = K_CT; K = L_CT
%K = K_CT; L = L_CT
K_FF = eye(2);
%K_FF = [1000 1000; 10000 10000]* 1/100000
%K_FF = inv(C*inv(eye(5)-Ad+Bd*L)*Bd)
% b)

% Qx(6,6) = 1/(1000^2); 
% Qu %Qu = zeros(1,1); Qu = 1/(10^2)*eye(2); 
%% Extend
C = [0 1 0 0 1]
Ad
Ad_e = [Ad zeros(5,1); -C 1]
Bd_e = [Bd ; 0 0]
%Cd_e = [C(1,:) 0; C(2,:) 0];
Cd_e =[C(1,:) 0];
sysextended = ss(Ad_e,Bd_e,Cd_e, 0, h);
[Klqr, S,e] = dlqr(Ad_e, Bd_e, Qx, Qu); 
K_fb = Klqr(:,1:length(Ad))
%K_I = Klqr(:,length(Ad)+1:end) 
K_I = Klqr(:,6)
Klqr