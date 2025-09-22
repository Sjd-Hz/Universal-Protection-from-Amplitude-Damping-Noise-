% Clear workspace and declare symbolic variables
clear all; clc; close all;

syms alpha beta lambda1 lambda2 t real
assume([t lambda1 lambda2], 'real');
assume(0 <= t & t <= 1);
assume(0 <= lambda1 & lambda1 <= 1);
assume(0 <= lambda2 & lambda2 <= 1);
assume(alpha^2+beta^2==1)
% Define computational basis states
H = [1; 0];  % |0>
V = [0; 1];  % |1>

% Initial Bell state |Phi+> = (|00> + |11>)/sqrt(2)
psi_ab = (1/sqrt(2)) * (kron(H, H) + kron(V, V));
rho_c = psi_ab * psi_ab';  % Density matrix of the Bell state

% Define amplitude damping coefficients
r1 = 1 - exp(-lambda1 * t);
r2 = 1 - exp(-lambda2 * t);

% Define ADC Kraus operators for each qubit
E0_r1 = [1, 0; 0, sqrt(1 - r1)];
E1_r1 = [0, sqrt(r1); 0, 0];
K1 = {E0_r1, E1_r1};

E0_r2 = [1, 0; 0, sqrt(1 - r2)];
E1_r2 = [0, sqrt(r2); 0, 0];
K2 = {E0_r2, E1_r2};

% Apply ADC channel to both qubits of the Bell state
rho_bell_adc = zeros(size(rho_c));
for i = 1:2
    for j = 1:2
        K = kron(K1{i}, K2{j});
        rho_bell_adc = simplify(rho_bell_adc + K * rho_c * K');
    end
end
disp('Density matrix after ADC:');
disp(simplify(rho_bell_adc));

% Input state rho_in = |psi><psi| = (alpha|0> + beta|1>)(...)
rho_in = [alpha^2, alpha * beta; alpha * beta, beta^2];

% Combine input state with noisy Bell channel
rho_com = kron(rho_in, rho_bell_adc);

% Define Bell basis projectors for Alice's measurement
b1 = (1/sqrt(2)) * (kron(H, H) + kron(V, V));
B1 = kron(b1 * b1', eye(2));

b2 = (1/sqrt(2)) * (kron(H, H) - kron(V, V));
B2 = kron(b2 * b2', eye(2));

b3 = (1/sqrt(2)) * (kron(H, V) + kron(V, H));
B3 = kron(b3 * b3', eye(2));

b4 = (1/sqrt(2)) * (kron(H, V) - kron(V, H));
B4 = kron(b4 * b4', eye(2));

% Compute probability of successful outcomes
P1 = simplify(trace(B1 * rho_com));
P3 = simplify(trace(B3 * rho_com));

% Define Pauli matrices
sigmax = [0, 1; 1, 0];
sigmaz = [1, 0; 0, -1];

% Compute Bob's conditional states for each Bell measurement
rho_B1 = simplify(PartialTrace(B1 * rho_com * B1', [1, 2]));
rho_B2 = simplify(PartialTrace(B2 * rho_com * B2', [1, 2]));
rho_B3 = simplify(PartialTrace(B3 * rho_com * B3', [1, 2]));
rho_B4 = simplify(PartialTrace(B4 * rho_com * B4', [1, 2]));

% Apply correction operations on Bob's qubit
rho_U1 = rho_B1;
rho_U2 = sigmaz * rho_B2 * sigmaz';
rho_U3 = sigmax * rho_B3 * sigmax';
rho_U4 = sigmax * sigmaz * rho_B4 * (sigmax * sigmaz)';

% Normalize the corrected states
g1 = simplify(trace(rho_U1));

ho_U1_nor = simplify(rho_U1 / g1);

g2 = simplify(trace(rho_U2));

ho_U2_nor = simplify(rho_U2 / g2);

g3 = simplify(trace(rho_U3));

ho_U3_nor = simplify(rho_U3 / g3);

g4 = simplify(trace(rho_U4));

ho_U4_nor = simplify(rho_U4 / g4);

% Compute fidelities with original input state
fid1 = simplify(trace(rho_in * rho_U1));
fid3 = simplify(trace(rho_in * rho_U3));

% Compute total and average fidelity (symmetry used: P1=P2, P3=P4)
Fidtot = simplify(2 * fid1 + 2 * fid3);
Fid_av = simplify(int(Fidtot, beta, 0, 1));
disp('Average teleportation fidelity:');
disp(simplify(Fid_av));