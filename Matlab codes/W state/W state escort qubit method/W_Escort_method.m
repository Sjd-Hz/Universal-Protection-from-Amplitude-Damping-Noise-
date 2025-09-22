%% Symbolic Definitions
clc; clear;
syms r q al be l1 l2 l3 t r1 r2 r3 real
assume([r, q, al, be, l1, l2, l3, r1, r2, r3] >= 0 & ...
       [r, q, al, be, l1, l2, l3, r1, r2, r3] <= 1);
assume(al^2 + be^2 == 1);

%% Basic Qubit States and Operators
H = [1; 0];  % |0⟩
V = [0; 1];  % |1⟩
sigmax = [0 1; 1 0];
sigmaz = [1 0; 0 -1];
Id2 = eye(2);

%% Construct Original W State
psi_W = 1/2 * (kron(V, kron(H, H)) + ...
               kron(H, kron(V, H)) + ...
               sqrt(2) * kron(H, kron(H, V)));
rho_W = psi_W * psi_W';

%% Add Escort Qubits and Construct Total Initial State
escort = kron(V, kron(V, V));  % |111⟩
initial_state = kron(psi_W, escort);  % q1 q2 q3 e1 e2 e3
rho_w0 = initial_state * initial_state';

%% Apply Pre-CNOT (q_i → e_i)
U_pre = CNOT_on_nqubits(6, 1, 4) * ...
         CNOT_on_nqubits(6, 2, 5) * ...
         CNOT_on_nqubits(6, 3, 6);
rho_c = U_pre * rho_w0 * U_pre';

%% Define Kraus Operators for Amplitude Damping
r1 = 1 - exp(-l1 * t);
r2 = 1 - exp(-l2 * t);
r3 = 1 - exp(-l3 * t);

E0_r1 = [1, 0; 0, sqrt(1 - r1)];  E1_r1 = [0, sqrt(r1); 0, 0];  K1 = {E0_r1, E1_r1};
E0_r2 = [1, 0; 0, sqrt(1 - r2)];  E1_r2 = [0, sqrt(r2); 0, 0];  K2 = {E0_r2, E1_r2};
E0_r3 = [1, 0; 0, sqrt(1 - r3)];  E1_r3 = [0, sqrt(r3); 0, 0];  K3 = {E0_r3, E1_r3};

%% Apply ADC to 6-Qubit System
rho_w_adc = zeros(size(rho_c));
for i1 = 1:2, for i2 = 1:2, for i3 = 1:2
for j1 = 1:2, for j2 = 1:2, for j3 = 1:2
    Kq1 = K1{i1}; Kq2 = K2{i2}; Kq3 = K3{i3};
    Ke4 = K1{j1}; Ke5 = K2{j2}; Ke6 = K3{j3};
    K_op = kron(Kq1, kron(Kq2, kron(Kq3, kron(Ke4, kron(Ke5, Ke6)))));
    rho_w_adc = simplify(rho_w_adc + K_op * rho_c * K_op');
end, end, end, end, end, end

%% Display Non-Zero Diagonal Elements
diagonal_elements = diag(rho_w_adc);
fprintf('Non-zero diagonal elements of the density matrix after ADCs:\n');
for i = 1:length(diagonal_elements)
    if ~isequal(diagonal_elements(i), sym(0))
        fprintf('rho(%d,%d) = %s\n', i, i, char(diagonal_elements(i)));
    end
end

%% Apply Post-CNOT and Project onto |111⟩ Escort Qubits
U_post = U_pre;
rho_post = U_post * rho_w_adc * U_post';
P_escort = projector_on_escort_111([4, 5, 6], 6);
rho_success = P_escort * rho_post * P_escort';
p11 = trace(rho_success);

% Trace out escort qubits
rho_shared = PartialTrace(rho_success, [4, 5, 6]);
disp('Shared non-normalized entangled state:');
disp(simplify(rho_shared));

%% Input Qubit and Combined System
rho_in = [al^2, al*be; al*be, be^2];
rho_com = kron(rho_in, rho_shared);

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
Fid_av = int(fid_tot, be, 0, 1)
