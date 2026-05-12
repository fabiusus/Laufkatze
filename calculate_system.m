%% 0. Symbolische Variablen initialisieren
clear; clc;

syms l m_1 m_s M k_r r g a real
syms phi x xdot xdotdot phidot phidotdot real 
syms z1 z2 z3 z4 u real
syms y real

%% 1. Gleichungen aufstellen und nach höchsten Ableitungen auflösen

R1 = phidotdot * l^2 * m_1 + xdotdot*l*m_1*cos(phi) + g*l*m_1*sin(phi) == 0;
R2 = l*m_1*sin(phi)*phidot^2 - xdotdot*m_1 - xdotdot*m_s - phidotdot*l*m_1*cos(phi) - k_r*xdot + M/r == 0;

% MATLAB löst das Gleichungssystem nach phi_ddot und x_ddot auf
sol = solve([R1, R2], [phidotdot, xdotdot]);

phidotdot_sol = simplify(sol.phidotdot);
xdotdot_sol   = simplify(sol.xdotdot);

%% 2. Nichtlineare Zustandsraumdarstellung aufstellen
% Wir substituieren die alten physikalischen Variablen durch die Zustandsvariablen
% z = [z1; z2; z3; z4] = [phi; phidot; x; xdot]
sub_old = [phi, phidot, x, xdot, M];
sub_new = [z1,  z2,     z3, z4,   u];

f_phidotdot = subs(phidotdot_sol, sub_old, sub_new);
f_xdotdot   = subs(xdotdot_sol,   sub_old, sub_new);

% Nichtlinearer Zustandsvektor dz/dt = f(z,u)
dz = [z2;            
      f_phidotdot;   
      z4;            
      f_xdotdot];    

dg = [z1;z3];

%% 3. Symbolische Linearisierung um die untere Ruhelage
% Pendel hängt nach unten (phi=0), alles steht still.
z_eq = [0; 0; 0; 0];
u_eq = 0;

% Jacobi-Matrix berechnen für A-Matrix (df/dz) und B-Matrix (df/du)
A_sym = jacobian(dz, [z1, z2, z3, z4]);
B_sym = jacobian(dz, u);
C_sym = jacobian(dg, [z1, z2, z3, z4]);
D_sym = jacobian(dg, u);

disp(A_sym);
disp(B_sym);

% Ruhelage in die Matrizen einsetzen
A_lin_sym = simplify(subs(A_sym, [z1, z2, z3, z4, u], [z_eq', u_eq]));
B_lin_sym = simplify(subs(B_sym, [z1, z2, z3, z4, u], [z_eq', u_eq]));
C_lin_sym = simplify(subs(C_sym, [z1, z2, z3, z4, u], [z_eq', u_eq]));

disp(A_lin_sym),
disp(B_lin_sym);
disp(C_lin_sym);

%% 4. Numerische Auswertung (Parameter einsetzen)
% Jetzt definieren wir deine Zahlenwerte. 
% Hinweis: m_l aus deinem Code wurde hier zu m_1 gemacht, damit es zur Formel passt.
vars_sym = [g,    m_s, m_1, l, k_r,    r];
vars_num = [9.81, 0.8, 0.1, 1, 7.2e-3, 0.041];

% Werte in die symbolischen A- und B-Matrizen einsetzen und in double umwandeln
A_num = double(subs(A_lin_sym, vars_sym, vars_num));
B_num = double(subs(B_lin_sym, vars_sym, vars_num));

disp('--- Numerische lineare Systemmatrix A ---');
disp(A_num);
disp('--- Numerische lineare Eingangsmatrix B ---');
disp(B_num);

%% 5. Stabilitätsanalyse über Eigenwerte
% Eigenwerte der numerischen A-Matrix berechnen
eig_A = eig(A_num);

disp('--- Eigenwerte der A-Matrix ---');
disp(eig_A);

% Automatische Auswertung der Stabilität anhand der Realteile der Eigenwerte
disp('--- Stabilitätsaussage ---');
% Wir nutzen eine kleine Toleranz (1e-6), um numerische Ungenauigkeiten abzufangen
max_real = max(real(eig_A));

if max_real < -1e-6
    disp('Ergebnis: Das System ist ASYMPTOTISCH STABIL (Alle Realteile < 0).');
elseif max_real > 1e-6
    disp('Ergebnis: Das System ist INSTABIL (Mindestens ein Realteil > 0).');
else
    disp('Ergebnis: Das System ist GRENZSTABIL (Größter Realteil ist 0).');
    disp('Grund: Der Kranwagen kann frei auf der x-Achse rollen (Position x ist unbestimmt).');
end

%% Steuerbarkeit und Beobachtbarkeit 

% Überprüfen der Steuerbarkeit
controllabilityMatrix = [B_num, A_num * B_num, A_num^2 * B_num, A_num^3 * B_num];
rankControllability = rank(controllabilityMatrix);

% Überprüfen der Steuerbarkeit
if rankControllability == size(A_num, 1)
    disp('Das System ist steuerbar.');
else
    disp('Das System ist nicht steuerbar.');
end

% Überprüfen der Beobachtbarkeit
observabilityMatrix = [C_sym; C_sym * A_sym; C_sym * A_sym^2; C_sym * A_sym^3];
rankObservability = rank(observabilityMatrix);

% Überprüfen der Beobachtbarkeit
if rankObservability == size(A_sym, 1)
    disp('Das System ist beobachtbar.');
else
    disp('Das System ist nicht beobachtbar.');
end


%% Übertragungsfunktionen
s = tf('s');

G_u_y1 = simplify(C_lin_sym(1,:) * inv(s * eye(size(A_num)) - A_num) * B_num)
G_u_y2 = simplify(C_lin_sym(2,:) * inv(s * eye(size(A_num)) - A_num) * B_num)


%% 4.2.3

function u = calculate_u(M, xdot)
    % Parameter
    P_max = 3.56; % [cite: 147]
    r = 0.041; % [cite: 170]
    
    % 1. Einzelmoment und omega berechnen
    M_M = M / 2; % 
    omega_M = -xdot / r; % 
    
    % Sonderfall abfangen (Division durch Null verhindern)
    if omega_M == 0
        omega_M = 1e-6; 
    end
    
    % 2. Vorzeichen und Beträge
    s = sign(M_M);
    if s == 0
        s = 1; % Verhindert Probleme bei M = 0
    end
    M_abs = abs(M_M);
    omega_tilde = s * omega_M;
    
    % 3. Fall A (Sättigung) testen
    P_sat = (9.434 * M_abs)^2; % 
    M_calc_test = (P_sat / abs(omega_M)) - sign(omega_tilde) * M_abs;
    
    if M_calc_test >= M_abs
        P_abs = P_sat;
    else
        % 4. Fall B (Linearer Bereich)
        x = (abs(omega_M) / 2) * ( (sign(omega_tilde)/9.434) + sqrt( (1/9.434^2) + (4*M_abs/abs(omega_M)) ) );
        P_abs = x^2;
    end
    
    % 5. Leistung berechnen und in Stellgröße umwandeln
    P = s * P_abs;
    u_raw = (100 * P) / P_max; % 
    
    % 6. Sättigung der Stellgröße auf +/- 100% [cite: 131, 134]
    u = max(-100, min(100, u_raw));
end

%% 4.2.4
u = 625 *sign(M)*M^2;

%% 7.ruhelagen
eq_ruhe = subs(f_phidotdot, {z2, z4, u}, {0, 0, 0}) == 0;
z1_sym = solve(eq_ruhe, z1); 
z1_num = double(subs(z1_sym(1), vars_sym, vars_num));
z_eq_num = [z1_num; 0; 0; 0];
disp('--- Ruhelagen der Zustandsvariablen (z_eq_num) ---');


