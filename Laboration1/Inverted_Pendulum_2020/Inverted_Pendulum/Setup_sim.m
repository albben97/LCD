
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
A = zeros(4,4);
% B
B = zeros(4,1);
% C
C = zeros(2,4);
% D
D = zeros(2,1);
%
%% Augment State-Space Matrices
% Ai
Ai = zeros(5,5);
% Bi
Bi = zeros(5,1);
% Ci
Ci = zeros(3,5);
% Di
Di = zeros(3,1);
%

%% DT State-Space
sysc = ss(A,B,C,D);
sysi = ss(Ai,Bi,Ci,Di);

%% Calculate LQR Control Gain
% Set LQR weighting matrices.
% Qx
Qx = zeros(5,5);
% Qu
Qu = zeros(1,1);
% Calculate control gain using LQR optimization.
% K_CT
K_CT = lqr(Ai,Bi,Qx,Qu);

%
%% Calculate Kalman Filter Gain
% Set Kalman Filter weighting matrices
% Qw
Qw = zeros(2,2);
% Qv
Qv = zeros(2,2);
% Nn
Nn = zeros(2,2);
% Calculate filter gain
% L_CT
G = zeros(4,2);
plant=ss(sysc.a,[sysc.b G],sysc.c,[sysc.d zeros(2)]);
% Calculate filter gain
[~,L_CT,~] = kalman(plant,Qw,Qv,Nn);
%
%% Display
disp( ' ' );
disp( 'SRV02+SIP Setup Succeed!' );
%
%%
