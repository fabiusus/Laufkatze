%% Initialisierung des Reglers
M_max = 0.4;
E_final = zeros(4,2);
Kx_final = zeros(1,4);
Kxi_final = zeros(1,1);

m1 = 1;
g=9.81;

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
K3place = place(A_erw, [B ; 0], pole);


Q = [diag([3,0.01,15,4,3]), diag([3,0.01,15,4,3]), diag([3,0.01,15,4,3]), diag([3,0.01,15,4,3]), diag([3,0.01,15,4,3]), diag([3,0.01,15,4,3])];
S = [1, 1, 1, 1, 1, 1];
Kx=zeros(6, 4);
Kxi=zeros(6,1);
for k = 1:6
 K_temp = lqr(A_erw,[B;0],Q(:,(((k-1)*5+1):k*5)),S(k));
 Kx(k,:) = K_temp(:,1:4);
 Kxi(k,:) = -K_temp(:,5);
end

zhut0 = [0;0;0;0];
xi0 = 0;

% E = place(A',C_new', 2*pole_ne)';

poles = [eig(A-B*Kx(1,:)), eig(A-B*Kx(2,:)), eig(A - B*Kx(3,:)), eig(A - B*Kx(4,:)), eig(A - B*Kx(5,:)), eig(A - B*Kx(6,:))];
E_temp = zeros(4,12);
for i = 1:6
    E_temp(:,2*i-1:2*i) = place(A', C_new', 2*poles(:,i))';
end

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


% 
% %%jold stuff
% 
% Q1 = diag([3,0.01,15,4,3]);
% S1 = diag(1);
% K31 = lqr(A_erw,[B;0],Q,S);
% Kx1 = K3(:,1:4);
% Kxi1 = -K3(:,5);
% 
% Q2 = diag([3,0.01,15,4,3]);
% S2 = diag(1);
% K32 = lqr(A_erw,[B;0],Q,S);
% Kx2 = K3(:,1:4);
% Kxi2 = -K3(:,5);
% 
% Q3 = diag([3,0.01,15,4,3]);
% S3 = diag(1);
% K33 = lqr(A_erw,[B;0],Q,S);
% Kx3 = K3(:,1:4);
% Kxi3 = -K3(:,5);
% 
% Q4 = diag([3,0.01,15,4,3]);
% S4 = diag(1);
% K34 = lqr(A_erw,[B;0],Q,S);
% Kx4 = K3(:,1:4);
% Kxi4 = -K3(:,5);
% 
% Q5 = diag([3,0.01,15,4,3]);
% S5 = diag(1);
% K35 = lqr(A_erw,[B;0],Q,S);
% Kx5 = K3(:,1:4);
% Kxi5 = -K3(:,5);
% 
% Q6 = diag([3,0.01,15,4,3]);
% S6 = diag(1);
% K36 = lqr(A_erw,[B;0],Q,S);
% Kx6 = K3(:,1:4);
% Kxi6 = -K3(:,5);