clc; clear;

% Define symbolic parameters
syms q al be l1 l2 l3 l t q1 q2 q3 q  v1 v2 v3 v real
assume([ q, al, be, l1, l2, l3, l, t, q1, q2, q3, q, v1, v2, v3, v] >= 0 & ...
       [ q, al, be, l1, l2, l3, l, t, q1, q2, q3, q, v1, v2, v3, v] <= 1)
assume(al^2 + be^2 == 1);

% Pauli matrices and Identity
sigmax = [0 1; 1 0];
sigmaz = [1 0; 0 -1];
Id2 = eye(2); % 2x2 Identity

% Choose mode: true = equal params, false = non-equal
equal_params = true;

% W-state definition
H = [1; 0]; V = [0; 1];
psi_W = 1/2 * (kron(V, kron(H, H)) + kron(H, kron(V, H)) + sqrt(2) * kron(H, kron(H, V)));
rho_W = psi_W * psi_W';

%% Weak Measurement (WM)
if equal_params
    W = [1, 0; 0, sqrt(1 - q1)];
    WM = kron(W, kron(W, W));
else
    W1 = [1, 0; 0, sqrt(1 - q1)];
    W2 = [1, 0; 0, sqrt(1 - q2)];
    W3 = [1, 0; 0, sqrt(1 - q3)];
    WM = kron(W1, kron(W2, W3));
end
rho_WM = simplify(WM * rho_W * WM');
disp('Density matrix after weak measurement operations:');
disp(simplify(rho_WM));
%% Amplitude Damping Channel (ADC)
if equal_params
    r1 = 1 - exp(-l1 * t);
    r2 = r1;
    r3 = r1;
else
    r1 = 1 - exp(-l1 * t);
    r2 = 1 - exp(-l2 * t);
    r3 = 1 - exp(-l3 * t);
end

K1 = {[1, 0; 0, sqrt(1 - r1)], [0, sqrt(r1); 0, 0]};
K2 = {[1, 0; 0, sqrt(1 - r2)], [0, sqrt(r2); 0, 0]};
K3 = {[1, 0; 0, sqrt(1 - r3)], [0, sqrt(r3); 0, 0]};

rho_ADC = zeros(size(rho_W));
for i = 1:2
    for j = 1:2
        for k = 1:2
            K = kron(K1{i}, kron(K2{j}, K3{k}));
            rho_ADC = simplify(rho_ADC + K * rho_WM * K');
        end
    end
end
disp('Density matrix after ADCs:');
disp(simplify(rho_ADC));
%% Reversal Weak Measurement (RWM)
if equal_params
    v1 = simplify(q1 + r1 * (1 - q1));
    M = [sqrt(1 - v1), 0; 0, 1];
    MR = kron(M, kron(M, M));
else
    v1 = simplify(q1 + r1 * (1 - q1));
    v2 = simplify(q2 + r2 * (1 - q2));
    v3 = simplify(q3 + r3 * (1 - q3));
    M1 = [sqrt(1 - v1), 0; 0, 1];
    M2 = [sqrt(1 - v2), 0; 0, 1];
    M3 = [sqrt(1 - v3), 0; 0, 1];
    MR = kron(M1, kron(M2, M3));
end

rho_RWM = simplify(MR * rho_ADC * MR');
prob = simplify(trace(rho_RWM))  % Overall success probability
disp('Density matrix after measurement reversal:');
disp(simplify(rho_RWM));

%% Teleportation process begins here %%%

% Define unknown input state to teleport
rho_in = [al^2, al*be; al*be, be^2];

% Combine input and shared entangled state
rho_com = kron(rho_in, rho_RWM);

% Define joint measurement projectors (Alice’s side for W state)
psi_P1 = 1/2 * (kron(H, kron(V, H)) + kron(H, kron(H, V)) + sqrt(2)*kron(V, kron(H, H)));
psi_P2 = 1/2 * (kron(H, kron(V, H)) + kron(H, kron(H, V)) - sqrt(2)*kron(V, kron(H, H)));
psi_P3 = 1/2 * (kron(V, kron(V, H)) + kron(V, kron(H, V)) + sqrt(2)*kron(H, kron(H, H)));
psi_P4 = 1/2 * (kron(V, kron(V, H)) + kron(V, kron(H, V)) - sqrt(2)*kron(H, kron(H, H)));
Proj1 = kron(psi_P1 * psi_P1', Id2);
Proj2 = kron(psi_P2 * psi_P2', Id2);
Proj3 = kron(psi_P3 * psi_P3', Id2);
Proj4 = kron(psi_P4 * psi_P4', Id2);

% Bob's state after each measurement (partial trace over Alice’s qubits)
rho_Proj1 = simplify(PartialTrace(Proj1 * rho_com * Proj1', [1, 2, 3]));
rho_Proj2 = simplify(PartialTrace(Proj2 * rho_com * Proj2', [1, 2, 3]));
rho_Proj3 = simplify(PartialTrace(Proj3 * rho_com * Proj3', [1, 2, 3]));
rho_Proj4 = simplify(PartialTrace(Proj4 * rho_com * Proj4', [1, 2, 3]));

% Bob’s recovery unitaries
U1 = Id2;
U2 = sigmaz;
U3 = sigmax;
U4 = sigmax * sigmaz;

% Apply recovery unitaries
rho_U1 = U1 * rho_Proj1 * U1';
rho_U2 = U2 * rho_Proj2 * U2';
rho_U3 = U3 * rho_Proj3 * U3';
rho_U4 = U4 * rho_Proj4 * U4';

% Probabilities of each outcome
g1 = simplify(trace(rho_U1));
g2 = trace(rho_U2);
g3 = trace(rho_U3);
g4 = trace(rho_U4);

% Total teleportation success probability (only for outcomes 1 and 3)
g_tot = simplify(2 * (g1 + g3))

% Normalize the output states
rho_U1_nor = rho_U1 / g1;
rho_U2_nor = rho_U2 / g2;
rho_U3_nor = rho_U3 / g3;
rho_U4_nor = rho_U4 / g4;

% Fidelity calculations between input and output states
fid1 = simplify(trace(rho_in * rho_U1_nor));
fid2 = simplify(trace(rho_in * rho_U2_nor));
fid3 = simplify(trace(rho_in * rho_U3_nor));
fid4 = simplify(trace(rho_in * rho_U4_nor));

% Average teleportation fidelity (based on success outcomes)
fid_tot = simplify(2 * (g1 * fid1 + g3 * fid3));

% Integrate over β to average over all input states on the Bloch sphere
Fid_av = int(fid_tot, be, 0, 1)
