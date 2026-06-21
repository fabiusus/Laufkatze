%% Initialisierung des Reglers
clear all
clc
% tunable parameters 

T = [3,3,1.5,2,2.3,3];
m1_sim = 1.1;
x_soll = 0.2;
Q = [diag([5000,1,500,400,0.1,10,1,10]), diag([5000,1,500,400,0.1,10,1,10]), diag([10000,0.01,5,800,0.01,100,0.01,100]), diag([10000,0.001,50,800,0.01,100,0.01,100]), diag([10000,0.01,5,800,0.01,100,0.01,100]), diag([10000,0.01,5,800,0.01,100,0.01,100])];
S = [3000, 3000, 5000, 20, 20, 20];
obs_mult = [2, 2, 5.0, 5.0, 7.0, 7.0];

G = [1;1;1;1];
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

B_erww = [B;0;0;0;0];

Kx=zeros(6, 4);
Kxi=zeros(6,4);
E_temp = zeros(4, 12);


zhut0 = [0;0;0;0];
xi0 = [0;0;0;0];
Q_sy = [0 1 0 0];
Pi_temp = zeros(4, 24);   
Gamma_temp = zeros(6, 4); 
for k = 1:6
    T_k = T(k);
    m1_k = m1_vec(k);
    A_k = [0 1 0 0;
        -g*(1+m1_k/0.8) 0 0 0.009;
        0 0 0 1;
        12.26*m1_k 0 0 -0.009];
    S_sy=[0 0 0 0;
        0 0 1 0;
        0 0 0 1;
        1/T_k^3 -1/T_k^3 -3/T_k^2 -3/T_k];
    A_erw = [A_k zeros(4,4);
            -G*C_vorfilter S_sy];
            
   
    Q_k = Q(:, (((k-1)*8+1):k*8));  
    K_temp = lqr(A_erw, B_erww, Q_k, S(k));
    Kx(k,:) = -K_temp(:,1:4);
    Kxi(k,:) = -K_temp(:,5:8);

    
    
    poles_k = eig(A_k + B*Kx(k,:));
    E_temp(:, 2*k-1:2*k) = place(A_k', C_new', obs_mult(k) * poles_k)';

    ff=[kron(S_sy',eye(4))-kron(eye(4),A_k), -kron(eye(4),B); kron(eye(4),C_vorfilter), zeros(4,4)]\[zeros(16,1);Q_sy(:)];
   
    Pi_temp(:, 4*k-3:4*k) = reshape(ff(1:16), 4, 4);
    Gamma_temp(k, :) = reshape(ff(17:end), 1, 4);
end

m1 = m1_sim;
A = [0 1 0 0;
    -g*(1+m1/0.8) 0 0 0.009;
    0 0 0 1;
    12.26*m1 0 0 -0.009];

