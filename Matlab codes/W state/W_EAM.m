%% Initialization and Symbolic Declarations
clc; clear;

syms r q alpha beta lambda1 lambda2 lambda3 t r1 r2 real
assume([r, q, alpha, beta, lambda1, lambda2, lambda3, r1, r2] >= 0 & ...
       [r, q, alpha, beta, lambda1, lambda2, lambda3, r1, r2] <= 1);
assume(alpha^2 + beta^2 == 1);

%% Define Qubit States and Operators
ket0 = [1; 0];  % |0⟩
ket1 = [0; 1];  % |1⟩
Id2 = eye(2);   % Identity matrix

%% Define W Entangled State
psi_W = 1/2 * (kron(ket1, kron(ket0, ket0)) + ...
               kron(ket0, kron(ket1, ket0)) + ...
               sqrt(2) * kron(ket0, kron(ket0, ket1)));
rho_W = psi_W * psi_W';

%% ADC Parameters (Decay Probabilities)
r1 = 1 - exp(-lambda1 * t);
r2 = 1 - exp(-lambda2 * t);
r3 = 1 - exp(-lambda3 * t);

%% Define ADC Kraus Operators (No-Jump Only)
E0_r1 = [1, 0; 0, sqrt(1 - r1)];
E0_r2 = [1, 0; 0, sqrt(1 - r2)];
E0_r3 = [1, 0; 0, sqrt(1 - r3)];

Kraus_no_jump = kron(E0_r1, kron(E0_r2, E0_r3));

%% Apply ADC to W State
rho_w_adc = Kraus_no_jump * rho_W * Kraus_no_jump';

%% Define Measurement Reversal Operators
v1 = r1;
v2 = r2;
v3 = r3;

M1 = [sqrt(1 - v1), 0; 0, 1];
M2 = [sqrt(1 - v2), 0; 0, 1];
M3 = [sqrt(1 - v3), 0; 0, 1];
M_total = kron(M1, kron(M2, M3));

%% Apply Reversal to the ADC Output
rho_after_reversal = simplify(M_total * rho_w_adc * M_total');
disp('Shared non-normalized entangled state:');
disp(simplify(rho_after_reversal));
%% Compute Probability and Normalize
prob = simplify(trace(rho_after_reversal));
rho_normalized = simplify(rho_after_reversal / prob);


%% Input Qubit and Combined System
rho_in = [alpha^2, alpha*beta; alpha*beta, beta^2];
rho_com = kron(rho_in, rho_after_reversal);

%% 
H = [1; 0];  % |0⟩
V = [0; 1];  % |1⟩
sigmax = [0 1; 1 0];
sigmaz = [1 0; 0 -1];
%% Define Alice's Measurement (W Basis)

psi_P1 = 1/2 * (kron(H, kron(V, H)) + kron(H, kron(H, V)) + sqrt(2)*kron(V, kron(H, H)));
psi_P2 = 1/2 * (kron(H, kron(V, H)) + kron(H, kron(H, V)) - sqrt(2)*kron(V, kron(H, H)));
psi_P3 = 1/2 * (kron(V, kron(V, H)) + kron(V, kron(H, V)) + sqrt(2)*kron(H, kron(H, H)));
psi_P4 = 1/2 * (kron(V, kron(V, H)) + kron(V, kron(H, V)) - sqrt(2)*kron(H, kron(H, H)));

rho_P1 = psi_P1 * psi_P1';
rho_P2 = psi_P2 * psi_P2';
rho_P3 = psi_P3 * psi_P3';
rho_P4 = psi_P4 * psi_P4';

Proj1 = kron(rho_P1, Id2);
Proj2 = kron(rho_P2, Id2);
Proj3 = kron(rho_P3, Id2);
Proj4 = kron(rho_P4, Id2);

%% Probabilities and Bob's State for Each Outcome
P1 = simplify(trace(Proj1 * rho_com * Proj1' / trace(rho_com)));
P2 = simplify(trace(Proj2 * rho_com * Proj2' / trace(rho_com)));
P3 = simplify(trace(Proj3 * rho_com * Proj3' / trace(rho_com)));
P4 = simplify(trace(Proj4 * rho_com * Proj4' / trace(rho_com)));

rho_Proj1 = simplify(PartialTrace(Proj1 * rho_com * Proj1', [1, 2, 3]));
rho_Proj2 = simplify(PartialTrace(Proj2 * rho_com * Proj2', [1, 2, 3]));
rho_Proj3 = simplify(PartialTrace(Proj3 * rho_com * Proj3', [1, 2, 3]));
rho_Proj4 = simplify(PartialTrace(Proj4 * rho_com * Proj4', [1, 2, 3]));

%% Bob's Recovery Operations
U1 = Id2;
U2 = sigmaz;
U3 = sigmax;
U4 = sigmax * sigmaz;

rho_U1 = U1 * rho_Proj1 * U1'; g1 = trace(rho_U1);
rho_U2 = U2 * rho_Proj2 * U2'; g2 = trace(rho_U2);
rho_U3 = U3 * rho_Proj3 * U3'; g3 = trace(rho_U3);
rho_U4 = U4 * rho_Proj4 * U4'; g4 = trace(rho_U4);

%% Total Teleportation Success Probability
g_tot = simplify(2 * (g1 + g3)) 

%% Normalized Output States
rho_U1_nor = rho_U1 / g1;
rho_U2_nor = rho_U2 / g2;
rho_U3_nor = rho_U3 / g3;
rho_U4_nor = rho_U4 / g4;

%% Teleportation Fidelities
fid1 = simplify(trace(rho_in * rho_U1_nor));
fid2 = simplify(trace(rho_in * rho_U2_nor));
fid3 = simplify(trace(rho_in * rho_U3_nor));
fid4 = simplify(trace(rho_in * rho_U4_nor));

%% Average Teleportation Fidelity
fid_tot = simplify(2 * (P1 * fid1 + P3 * fid3))
Fid_av = int(fid_tot, beta, 0, 1)

