%% Initialisierung des Reglers
clear all
clc
% tunable parameters 

T = [3,4,5,6,5,8];
m1_sim = 1;
x_soll = 0.2;
Q = [diag([10,1,10,10,1]), diag([1,1,1,1,1]), diag([1,1,1,1,1]), diag([1,1,1,1,1]), diag([1,1,1,1,1]), diag([1,1,1,1,1])];
S = [100, 1, 1, 1, 1, 1];
obs_mult = [1.9, 3, 5.0, 5.0, 7.0, 7.0];
%[5000,1,500,400,0.1,10,1,10]

G = [1];
q_basis = [10000,0.01,5,800,0.01,100,0.01,100];

m1_vec = [1,1.5,2.5,3.5,4.5,5];

% fixed parameters

g=9.81;

E_final = zeros(4,2);
Kx_final = zeros(1,4);
Kxi_final = zeros(1,1);
% System
x0= [0; 0; 0; 0];

B=[0; -30.49; 0; 30.49];
C = eye(4);
C_new= [1 0 0 0;
    0 0 1 0];
C_vorfilter = [0 0 1 0];
D = zeros(4,1);


% A= [0 1 0 0;
%     -g*(1+m1/0.8) 0 0 0.009;
%     0 0 0 1;
%     12.26*m1 0 0 -0.009];
% A_erw = [A zeros(4,4);
%         -G*C_vorfilter S_sy];

B_erww = [B;0];

Kx=zeros(6, 4);
Kxi=zeros(6,4);
E_temp = zeros(4, 12);


zhut0 = [0;0;0;0];
xihut0 = zeros(8,1);
xi0 = [0;0;0;0];
Q_sy = [1];
Pi_temp = zeros(4, 6);   
Gamma_temp = zeros(1, 6); 
for k = 1:6
    T_k = T(k);
    m1_k = m1_vec(k);
    A_k = [0 1 0 0;
        -g*(1+m1_k/0.8) 0 0 0.009;
        0 0 0 1;
        12.26*m1_k 0 0 -0.009];
    S_sy=0;
    A_erw = [A_k zeros(4,1);
            -G*C_vorfilter S_sy];
   
    Q_k = Q(:, (((k-1)*5+1):k*5));  
    K_temp = lqr(A_erw, B_erww, Q_k, S(k));
    Kx(k,:) = -K_temp(:,1:4);
    Kxi(k,:) = -K_temp(:,5);

    poles_k = eig(A_k + B*Kx(k,:));
    % A_beo = [A_k, ones(4,4);
    %     zeros(4,4), S_sy];
    
    E_temp(:, 2*k-1:2*k) = place(A_k', C_new', obs_mult(k) * poles_k)';
    % E_temp(:,2*k-1:2*k) = place(A_beo', [C_new zeros(2,4)]', obs_mult(k)*[poles_k; -2; -2.1; -2.2; -2.3])';

    ff=[kron(S_sy',eye(4))-kron(eye(1),A_k), -kron(eye(1),B); kron(eye(1),C_vorfilter), 0]\[zeros(4,1);Q_sy(:)];
   
    Pi_temp(:, k) = ff(1:4,1);
    Gamma_temp(k) = ff(5,1);
end

m1 = m1_sim;
A = [0 1 0 0;
    -g*(1+m1/0.8) 0 0 0.009;
    0 0 0 1;
    12.26*m1 0 0 -0.009];

