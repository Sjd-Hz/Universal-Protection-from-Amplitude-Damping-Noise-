clc; clear all;

% Define symbolic variables
syms r q alpha beta lambda1 lambda2 lambda3 t r1 r2 r3 real
assume([r, q, alpha, beta, lambda1, lambda2, lambda3] >= 0 & ...
       [r, q, alpha, beta, lambda1, lambda2, lambda3] <= 1 & alpha^2 + beta^2 == 1)

%% Define initial basis and operators
H = [1; 0]; % |0>
V = [0; 1]; % |1>
sigmax = [0 1; 1 0]; % Pauli-X
sigmaz = [1 0; 0 -1]; % Pauli-Z
I2 = eye(2); % 2x2 Identity

% W state (three-qubit entangled state)
psi_W = 1/2 * (kron(V, kron(H, H)) + kron(H, kron(V, H)) + sqrt(2) * kron(H, kron(H, V)));
rho_W = psi_W * psi_W';

%% Amplitude Damping Channel (ADC)

% Define damping probabilities
r1 = 1 - exp(-lambda1 * t);
r2 = 1 - exp(-lambda2 * t);
r3 = 1 - exp(-lambda3 * t);

% Define ADC Kraus operators for each qubit
E0_1 = [1, 0; 0, sqrt(1 - r1)]; E1_1 = [0, sqrt(r1); 0, 0];
E0_2 = [1, 0; 0, sqrt(1 - r2)]; E1_2 = [0, sqrt(r2); 0, 0];
E0_3 = [1, 0; 0, sqrt(1 - r3)]; E1_3 = [0, sqrt(r3); 0, 0];

K1 = {E0_1, E1_1}; % ADC on qubit 1
K2 = {E0_2, E1_2}; % ADC on qubit 2
K3 = {E0_3, E1_3}; % ADC on qubit 3

% Apply ADC to W state
rho_W_ADC = zeros(size(rho_W));
for i = 1:2
    for j = 1:2
        for k = 1:2
            K = kron(K1{i}, kron(K2{j}, K3{k}));
            rho_W_ADC = simplify(rho_W_ADC + K * rho_W * K');
        end
    end
end

disp('Density matrix after ADC:');
disp(simplify(rho_W_ADC));

%% Input state and combined state
rho_in = [alpha^2, alpha*beta; alpha*beta, beta^2]; % pure input state
rho_combined = kron(rho_in, rho_W_ADC); % input ⊗ noisy W state

%% Define joint measurement basis (Alice's measurements)

psi_P1 = 1/2*(kron(H, kron(V, H)) + kron(H, kron(H, V)) + sqrt(2)*kron(V, kron(H, H)));
psi_P2 = 1/2*(kron(H, kron(V, H)) + kron(H, kron(H, V)) - sqrt(2)*kron(V, kron(H, H)));
psi_P3 = 1/2*(kron(V, kron(V, H)) + kron(V, kron(H, V)) + sqrt(2)*kron(H, kron(H, H)));
psi_P4 = 1/2*(kron(V, kron(V, H)) + kron(V, kron(H, V)) - sqrt(2)*kron(H, kron(H, H)));

rho_P1 = psi_P1 * psi_P1';
rho_P2 = psi_P2 * psi_P2';
rho_P3 = psi_P3 * psi_P3';
rho_P4 = psi_P4 * psi_P4';

Proj1 = kron(rho_P1, I2);
Proj2 = kron(rho_P2, I2);
Proj3 = kron(rho_P3, I2);
Proj4 = kron(rho_P4, I2);

P1 = simplify(trace(Proj1 * rho_combined * Proj1') / trace(rho_combined));
P2 = simplify(trace(Proj2 * rho_combined * Proj2') / trace(rho_combined));
P3 = simplify(trace(Proj3 * rho_combined * Proj3') / trace(rho_combined));
P4 = simplify(trace(Proj4 * rho_combined * Proj4') / trace(rho_combined));

%% Bob's conditional states
rho1 = simplify(PartialTrace(Proj1 * rho_combined * Proj1', [1,2,3]));
rho2 = simplify(PartialTrace(Proj2 * rho_combined * Proj2', [1,2,3]));
rho3 = simplify(PartialTrace(Proj3 * rho_combined * Proj3', [1,2,3]));
rho4 = simplify(PartialTrace(Proj4 * rho_combined * Proj4', [1,2,3]));

% Unitary corrections
U1 = I2;
U2 = sigmaz;
U3 = sigmax;
U4 = sigmax * sigmaz;

rho_out1 = U1 * rho1 * U1';
rho_out2 = U2 * rho2 * U2';
rho_out3 = U3 * rho3 * U3';
rho_out4 = U4 * rho4 * U4';

% Success probabilities
g1 = simplify(trace(rho_out1));
g2 = simplify(trace(rho_out2));
g3 = simplify(trace(rho_out3));
g4 = simplify(trace(rho_out4));

g_total = simplify(2 * (g1 + g3));
disp('Total teleportation success probability:');
disp(simplify(g_total));

% Normalized output states
rho1_nor = rho_out1 / g1;
rho2_nor = rho_out2 / g2;
rho3_nor = rho_out3 / g3;
rho4_nor = rho_out4 / g4;

% Fidelities
fid1 = simplify(trace(rho_in * rho1_nor));
fid2 = simplify(trace(rho_in * rho2_nor));
fid3 = simplify(trace(rho_in * rho3_nor));
fid4 = simplify(trace(rho_in * rho4_nor));

% Average teleportation fidelity
fid_total = simplify(2 * (P1 * fid1 + P3 * fid3));
Fid_av = simplify(int(fid_total, beta, 0, 1));
disp('Average teleportation fidelity:');
disp(simplify(Fid_av));

syms lambda
% Replace lambda1 = lambda2 = lambda3 = lambda
fid_total_equal_lambda = subs(fid_total, [lambda1, lambda2, lambda3], [lambda, lambda, lambda]);

% Replace alpha using alpha^2 + beta^2 = 1 → alpha = sqrt(1 - beta^2)
fid_total_equal_lambda = subs(fid_total_equal_lambda, alpha, sqrt(1 - beta^2));

% Recompute the average fidelity
Fid_av_equal_lambda = simplify(int(fid_total_equal_lambda, beta, 0, 1));

disp('Average teleportation fidelity with λ₁ = λ₂ = λ₃ = λ:');
disp(Fid_av_equal_lambda);