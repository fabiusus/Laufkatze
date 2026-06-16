clear all
close all
clc

m1 = 1;
g=9.81;
T=10
A= [0 1 0 0;
    -g*(1+m1/0.8) 0 0 0.009;
    0 0 0 1;
    12.26*m1 0 0 -0.009];
B=[0; -30.49; 0; 30.49];

C = eye(4);
C_new= [1 0 0 0;
    0 0 1 0];
C_vorfilter = [0 0 1 0];

D = zeros(4,1);

pole = -1*[2.5,2.6,2.3,2.4,10.1];
% pole = -1*[1.1,1.2,1.3,1.4];
pole_ne = -1*[2.5,2.6,2.3,2.4];
A_erw = [A zeros(4,1);
        -[0 0 1 0] zeros(1,1)];
x0= [0; 0; 1; 0];
% K = place(A,B,pole_ne);
K3place = place(A_erw, [B ; 0], pole)
Q= diag([3,0.01,15,4,3]);
S = diag(1);
K3 = lqr(A_erw,[B;0],Q,S)
Kx= K3(:,1:4);
Kxi = -K3(:,5);



zhut0 = [0;0;0;0];
xi0 = 0;
E = place(A',C_new', 2*pole_ne)';

% v = -1/(C_vorfilter*(inv(A-B*K)*B));
%%Sylvester Gleichingen
S_sy=[0 0 0 0;
    0 0 1 0;
    0 0 0 1;
    1/T^3 -1/T^3 -3/T^2 -3/T];
Q_sy = [0 1 0 0];
ff=[kron(S_sy',eye(4))-kron(eye(4),A), -kron(eye(4),B); -kron(eye(4),C_vorfilter), zeros(4,4)]\[zeros(16,1);Q_sy(:)];
pi_vec=ff(1:16);
gamma_vec=ff(17:end);
Pi=reshape(pi_vec,4,4);
Gamma=reshape(gamma_vec,1,4);
