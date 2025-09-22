clear all;
close all;
clc;

% Define symbolic parameters
syms q v lambda t al be real
assume(0 <= [q v lambda t al be] & [q v lambda  t al be] <= 1 & al^2 + be^2 == 1)

% Pauli operators
sigmax = [0 1; 1 0];
sigmaz = [1 0; 0 -1];
Id = eye(2);

% Define standard basis vectors
H = [1; 0];  % |0⟩
V = [0; 1];  % |1⟩

% Shared Bell state (before noise)
psi_bell = 1/sqrt(2) * (kron(H, H) + kron(V, V));
rho_bell = psi_bell * psi_bell';

% --- Weak measurement (WM) ---
W1 = [1, 0; 0, sqrt(1 - q)];
W2 = W1;
WM = kron(W1, W2);
rho_after_WM = simplify(WM * rho_bell * WM');

% --- Amplitude Damping Channel (ADC) ---
r = 1 - exp(-lambda * t);
E0 = [1, 0; 0, sqrt(1 - r)];
E1 = [0, sqrt(r); 0, 0];
K1 = {E0, E1};
K2 = {E0, E1};

rho_after_ADC = zeros(4, 4);
for i1 = 1:2
    for i2 = 1:2
        Kraus_op = kron(K1{i1}, K2{i2});
        rho_after_ADC = simplify(rho_after_ADC + Kraus_op * rho_after_WM * Kraus_op');
    end
end
disp('Density matrix after ADC:');
disp(simplify(rho_after_ADC));

% --- Reversal Weak Measurement (RWM) ---
v = q + r * (1 - q);  % optimized reversal strength
M1 = [sqrt(1 - v), 0; 0, 1];
M2 = M1;
RWM = kron(M1, M2);
rho_after_RWM = simplify(RWM * rho_after_ADC * RWM');
prob_success = simplify(trace(rho_after_RWM));  % overall success probability of WMR method

% Final shared state (non-normalized)
rho_shared = rho_after_RWM;

% --- Input state to be teleported ---
input_state = al * H + sqrt(1 - al^2) * V;
rho_in = simplify(input_state * input_state');

% Combine with shared entangled state
rho_combined = kron(rho_in, rho_shared);

% --- Bell state projectors (Alice's joint measurements) ---
b1 = 1/sqrt(2)*(kron(H, H) + kron(V, V));
b2 = 1/sqrt(2)*(kron(H, H) - kron(V, V));
b3 = 1/sqrt(2)*(kron(H, V) + kron(V, H));
b4 = 1/sqrt(2)*(kron(H, V) - kron(V, H));

B1 = kron(b1 * b1', Id);
B2 = kron(b2 * b2', Id);
B3 = kron(b3 * b3', Id);
B4 = kron(b4 * b4', Id);

% --- Bob's (unnormalized) states for each measurement outcome ---
rho_B1 = simplify(PartialTrace(B1 * rho_combined * B1', [1, 2]));
rho_B2 = simplify(PartialTrace(B2 * rho_combined * B2', [1, 2]));
rho_B3 = simplify(PartialTrace(B3 * rho_combined * B3', [1, 2]));
rho_B4 = simplify(PartialTrace(B4 * rho_combined * B4', [1, 2]));

% --- Apply correcting unitaries ---
U1 = Id;
U2 = sigmaz;
U3 = sigmax;
U4 = sigmax * sigmaz;

rho_U1 = simplify(U1 * rho_B1 * U1');
rho_U2 = simplify(U2 * rho_B2 * U2');
rho_U3 = simplify(U3 * rho_B3 * U3');
rho_U4 = simplify(U4 * rho_B4 * U4');

% --- Success probabilities ---
g1 = trace(rho_U1);
g2 = trace(rho_U2);
g3 = trace(rho_U3);
g4 = trace(rho_U4);
g_tot = simplify(2 * (g1 + g3));  % total useful teleportation probability
disp('Total teleportation success probability:');
disp(simplify(g_tot));
% --- Normalized states at Bob's side ---
rho_U1_nor = rho_U1 / g1;
rho_U2_nor = rho_U2 / g2;
rho_U3_nor = rho_U3 / g3;
rho_U4_nor = rho_U4 / g4;

% --- Fidelity calculations ---
fid1 = simplify(trace(rho_in * rho_U1_nor));
fid2 = simplify(trace(rho_in * rho_U2_nor));
fid3 = simplify(trace(rho_in * rho_U3_nor));
fid4 = simplify(trace(rho_in * rho_U4_nor));

% --- Average teleportation fidelity ---
P1 = simplify(trace(B1 * rho_combined / trace(rho_combined) * B1'));
P3 = simplify(trace(B3 * rho_combined / trace(rho_combined) * B3'));
fid_tot = simplify(2 * (P1 * fid1 + P3 * fid3));
Fid_av = simplify(int(fid_tot, al, 0, 1));
disp('Average teleportation fidelity:');
disp(simplify(Fid_av));