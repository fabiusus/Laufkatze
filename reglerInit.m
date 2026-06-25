%% Initialisierung des Reglers


T = [0.34, 0.34, 0.34, 0.42, 0.48, 0.46];

Q = { [900,8,800, 50,0.05], ...
      [900,8,800, 90,0.05], ...
      [900,8,800,140,0.05], ...
      [900,8,800,140,0.05], ...
      [900,8,800, 90,0.05], ...
      [900,8,1200,140,0.05] };
S = [600, 600, 600, 600, 600, 1200];
obs_mult = [9, 9, 9, 9, 9, 9];  

% ---- one-sided rail-safety guard params (Controller03Synth) ----
x_aim     = 2.43;     % reference setback: aim a few mm inside the rail
x_wall    = 2.44;     % virtual-wall trigger position
kwall     = 300.0;    % quadratic virtual-wall gain
kwall_lin = 15.0;     % linear virtual-wall gain
kdamp     = 8.0;      % outward-velocity damping gain
mclip     = 0.4;      % torque clip [N*m]  (== H.MMAX)

m1_vec = [1, 1.5, 2.5, 3.5, 4.5, 5.5];

% fixed parameters
g=9.81;
ms = 0.8;

A_grid     = zeros(4, 4, 6);
B_grid     = zeros(4, 1, 6);
Kx_grid    = zeros(1, 4, 6);
Kxi_grid   = zeros(1, 6);     % SCALAR integral gain per bin
E_grid     = zeros(4, 2, 6);
Pi_grid    = zeros(4, 6);
Gamma_grid = zeros(1, 6);
Apt3_grid  = zeros(3, 3, 6);
bpt3_grid  = zeros(3, 6);

% System

B=[0; -1/(0.041*1*ms); 0; 1/(0.041*ms)];
B_erww      = [B; 0];
C_new= [1 0 0 0;
    0 0 1 0];
C_vorfilter = [0 0 1 0];
G_int = 1;
S_sy = 0;
x_0 = [0;0;0;0];
xi0 = 0;
zhut0 = [0;0;0;0];


for k = 1:6
    m1 = m1_vec(k);
    A_k = [ 0,             1, 0, 0;
           -g*(1+m1/ms),   0, 0, 0.009;
            0,             0, 0, 1;
            g*m1/ms,       0, 0, -0.009 ];
    A_grid(:, :, k) = A_k;
    B_grid(:, :, k) = B;

    % --- LQI: augmented system with single integrator on position error ---
    A_erw  = [A_k, zeros(4, 1);
              -G_int*C_vorfilter, 0];
    K = lqr(A_erw, B_erww, diag(Q{k}), S(k));   % 1x5
    Kx  = -K(1:4);          % 1x4
    Kxi = -K(5);            % scalar
    Kx_grid(:, :, k) = Kx;
    Kxi_grid(k)      = Kxi;

    % --- Luenberger observer: poles = obs_mult * eig(A_k + B*Kx) ---
    poles = eig(A_k + B*Kx);
    E = place(A_k', C_new', obs_mult(k)*poles)';   % 4x2
    E_grid(:, :, k) = E;

    % --- steady-state feedforward: [-A_k,-B; C_vorfilter,0]\[0;0;0;0;1] ---
    ff = [-A_k, -B; C_vorfilter, 0] \ [0; 0; 0; 0; 1];
    Pi_grid(:, k)  = ff(1:4);
    Gamma_grid(k)  = ff(5);

    % --- PT3 prefilter (triple pole at -1/T), DC gain 1 ---
    Tk = T(k);
    Apt3_grid(:, :, k) = [ 0,         1,         0;
                           0,         0,         1;
                          -1/Tk^3,   -3/Tk^2,   -3/Tk ];
    bpt3_grid(:, k)    = [0; 0; 1/Tk^3];
end