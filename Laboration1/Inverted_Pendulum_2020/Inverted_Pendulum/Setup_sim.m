
%% Initialization Settings(DONT CHANGE ANYTHING !!)
clear all;
EXT_GEAR_CONFIG = 'HIGH';
ENCODER_TYPE = 'E';
TACH_OPTION = 'YES';
LOAD_TYPE = 'NONE';
K_AMP = 1;
AMP_TYPE = 'VoltPAQ';
VMAX_DAC = 10;
ROTPEN_OPTION = '2DGANTRY-E';
PEND_TYPE = 'MEDIUM_12IN';
[ Rm, kt, km, Kg, eta_g, Beq, Jm, Jeq_noload, eta_m, K_POT, K_TACH, K_ENC, VMAX_AMP, IMAX_AMP ] = config_srv02( EXT_GEAR_CONFIG, ENCODER_TYPE, TACH_OPTION, AMP_TYPE, LOAD_TYPE );
[ Lb, Jarm, K_POT_2DP, K_ENC_2DP ] = config_2d_gantry( Jeq_noload );
%
%%%
%%%%
%%%%%
%%%%%% DO NOT CHANGE ANYTHING ABOVE THIS AREA !!
%%%%
%%%
%
%% System State-Space Matrices
% A
%A = zeros(4,4);
A = [0         0    1.0000         0;
     0         0         0    1.0000;
     0  0.1*128.2092   -7.7744         0;
     0  0.1*528.8277 -4.5648         0];
% 528.8277
% B
%B = zeros(4,1);
B = [0;
     0;
     14.4784;
     8.5012];
% C
%C = zeros(2,4);
C = [1 0 0 0;
     0 1 0 0];
% D
%D = zeros(2,1);
D = [0; 0];
%

%% Augment State-Space Matrices
% Ai
%Ai = zeros(5,5);
% Bi
%Bi = zeros(5,1);
% Ci
%Ci = zeros(3,5);
% Di
%Di = zeros(3,1);
%
Ai = [A zeros(4,1); C(1,:) zeros(1,1)];

Bi = [B; 0];

Ci = [C, zeros(2,1); 0 0 0 0 1];

Di = [D; 0];

%% DT State-Space
sysc = ss(A,B,C,D);
sysi = ss(Ai,Bi,Ci,Di);

%% Calculate LQR Control Gain
% Set LQR weighting matrices.
% Qx
%Qx = zeros(5,5);
%Qx = 1/(20^2)*eye(5); %Simulation
Qx = diag([500 500 1 1 1]); % Laboration
%Qx = diag([500 500 0 0 1]); % Laboration 1.4
% Qu
%Qu = zeros(1,1);
%Qu = 1/(10^2); % Simulation
%Qu = 0.2; % Laboration
Qu = 0.2; % Laboration 1.4
% Calculate control gain using LQR optimization.
% K_CT
K_CT = lqr(Ai,Bi,Qx,Qu);

%
%% Calculate Kalman Filter Gain
% Set Kalman Filter weighting matrices
% Qw
%Qw = zeros(2,2);
%Qw = [10 0; 0 10]; % Simulation
Qw = 50*diag([1, 1]); % Laboration
% Qv
%Qv = [1 0; 0 1]; % Simulation
Qv = 0.0000001*diag([1,1]); % Laboration
% Nn
Nn = zeros(2,2);
%Nn = Qw*Qv.';
% Calculate filter gain
% L_CT
%G = zeros(4,2);
G = [0 0 1 0; 0 0 0 1]';
plant=ss(sysc.a,[sysc.b G],sysc.c,[sysc.d zeros(2)]);
% Calculate filter gain
[~,L_CT,~] = kalman(plant,Qw,Qv,Nn);
%
%% Display
disp( ' ' );
disp( 'SRV02+SIP Setup Succeed!' );
%
%%
% Original initial conditions: pi/180*[-5, 2, 4, 0.3]
% Laboration init conds:
init_condi = pi/180*[5, 5, 0, 0];
% Parameters
Mp = 0.1270;
Lr = 0.1270;
Lp = 0.3111;
Jr = 0.0083;
Jp = 0.0012;
Dr = 0.0690;
Dp = 0.024*0;
Co = 0.1285;
g = 0.981;