%% Laboration 1
% x1 = delta_theta, x2 = delta_alpha, x3 = delta_theta_dot, x4 = delta_alpha_dot 
% Parameters
Mp = 0.1270;
Lr = 0.1270;
Lp = 0.3111;
Jr = 0.0083;
Jp = 0.0012;
Dr = 0.0690;
Co = 0.1285;
g = 0.981;

% Exercise 1, calculated by hand
den = (4*Jr+4*Mp*Lr^2)*Jp + Mp*Lp^2*Jr

df1_dx1 = 0;
df1_dx2 = 0;
df1_dx3 = 1;
df1_dx4 = 0;
df2_dx1 = 0;
df2_dx2 = 0;
df2_dx3 = 0;
df2_dx4 = 1;
df3_dx1 = 0;
df3_dx2 = -1/2*(20*Mp*Lp^2*Lr*g)/den;
df3_dx3 = (-1/2*(8*Dr*Jp)+1/2*(-2*Mp*Lp^2*Dr))/den;
df3_dx4 = 0;
df4_dx1 = 0;
df4_dx2 = (20*Jr*Mp*Lp*g + 20*Mp^2*Lr^2*Lp*g)/den;
df4_dx3 = -2*Mp*Lr*Lp*Dr/den;
df4_dx4 = 0;


df1_du = 0;
df2_du = 0;
df3_du = (4*Co*Jp + Mp*Lp^2*Co)/den;
df4_du = (2*Mp*Lr*Lp*Co)/den;

% Exercise 2
A_delta = [df1_dx1 df1_dx2 df1_dx3 df1_dx4;
           df2_dx1 df2_dx2 df2_dx3 df2_dx4;
           df3_dx1 df3_dx2 df3_dx3 df3_dx4;
           df4_dx1 df4_dx2 df4_dx3 df4_dx4];

B_delta = [df1_du; df2_du; df3_du; df4_du];

C_delta = [1 0 0 0;
           0 1 0 0];
       
D_delta = [0; 0];

% Exercise 3
% Extended models, integrating action
Ae = [A_delta zeros(4,1); -C_delta(1,:) zeros(1,1)];

Be = [B_delta; 0];

Ce = [C_delta, zeros(2,1)];

De = D_delta;

% Exercise 4
% Eigenvalues to check stability
eig_A_delta = eig(A_delta)
eig_Ae = eig(Ae)

sys_delta = ss(A_delta,B_delta,C_delta,D_delta);
zeros_delta = tzero(sys_delta)

syse = ss(Ae,Be,Ce,De);
zerose = tzero(syse)

% The inclusion of an extra integral state does not change the stability
% property since we have just added a state which generates a pole and a
% zero that cancels each other. The systems are not assymptotically stable.
% The extended system has one zero (which is canceled by a pole) and the
% original delta-system has no zeroes.

% Exercise 5
% Observability
obsv_delta = obsv(sys_delta)
rank_obsv_delta = rank(obsv_delta) % full rank = 4
obsve = obsv(syse)
rank_obsve = rank(obsve) % full rank = 5

% Only the delta-system is observable, the extended integral state is not
% observable.

% Exercise 6
% OPTIONAL

% Exercise 7
% Controllability
ctrb_delta = ctrb(sys_delta)
rank_ctrb_delta = rank(ctrb_delta) % full rank = 4
ctrbe = ctrb(syse)
rank_ctrbe = rank(ctrbe) % full rank = 5

% Both systems are controllable.

% Exercise 8
% OPTIONAL

% Exercise 9
% Qu and Qx decided by Byron's rule
Qu = 1/10;
Qx_delta = 1/20*eye(4);
Qxe = 1/20*eye(5);

%K_delta = 1/Qu*B_delta'*
%Ke = 1/Qu*Be'

[K_delta, S_delta, e_delta] = lqr(sys_delta,Qx_delta,Qu);
[Ke, Se, ee] = lqr(syse,Qxe,Qu);

% Gain K for the delta-system as well as the extended system
K_delta
Ke


