clc; clear;

% Declare symbolic variables
syms r q al be l1 l2 t r1 r2 real
assume((0 <= r) & (r <= 1))
assume((0 <= l1) & (l1 <= 1))
assume((0 <= l2) & (l2 <= 1))
assume((0 <= r1) & (r1 <= 1))
assume((0 <= r2) & (r2 <= 1))
assume((0 <= q) & (q <= 1))
assume((0 <= al) & (al <= 1))
assume((0 <= be) & (be <= 1))

% Define |0> and |1> basis vectors
H = [1; 0]; % |0⟩
V = [0; 1]; % |1⟩

% Define Bell state: (|00⟩ + |11⟩) / sqrt(2)
bell_state = 1/sqrt(2) * (kron(H, H) + kron(V, V));

% Define escort qubits (ancillas), initialized in state |1⟩
ancilla1 = [0; 1];
ancilla2 = [0; 1];

% Construct the full initial state: Bell ⊗ anc1 ⊗ anc2
initial_state = kron(bell_state, kron(ancilla1, ancilla2));
rho_bell0 = initial_state * initial_state'; % Density matrix

% Define CNOT gate acting on qubit 1 (control) and 3 (target)
C13 = [1 0 0 0 0 0 0 0;
       0 1 0 0 0 0 0 0;
       0 0 1 0 0 0 0 0;
       0 0 0 1 0 0 0 0;
       0 0 0 0 0 1 0 0;
       0 0 0 0 1 0 0 0;
       0 0 0 0 0 0 0 1;
       0 0 0 0 0 0 1 0];

% Pre-noise CNOT gates (entangled qubits to escort qubits)
CNOT_1 = kron(C13, eye(2));  % CNOT between qubits 1 and 3
CNOT_2 = kron(eye(2), C13);  % CNOT between qubits 2 and 4
rho_c = CNOT_2 * CNOT_1 * rho_bell0 * CNOT_1' * CNOT_2'; % Apply both CNOTs

% Define ADC parameters
r1 = 1 - exp(-l1 * t);
r2 = 1 - exp(-l2 * t);

% Define Kraus operators for ADC on each qubit
E0_r1 = [1, 0; 0, sqrt(1 - r1)];
E1_r1 = [0, sqrt(r1); 0, 0];
K1 = {E0_r1, E1_r1};  % ADC on qubit 1 and escort 1

E0_r2 = [1, 0; 0, sqrt(1 - r2)];
E1_r2 = [0, sqrt(r2); 0, 0];
K2 = {E0_r2, E1_r2};  % ADC on qubit 2 and escort 2

% Apply ADC to the system (4 qubits total)
rho_bell_adc = zeros(size(rho_c));
for i1 = 1:2
    for i2 = 1:2
        for j1 = 1:2
            for j2 = 1:2
                % Tensor structure: K1 ⊗ K2 ⊗ K1 ⊗ K2
                K = kron(K1{i1}, kron(K2{j1}, kron(K1{i2}, K2{j2})));
                rho_bell_adc = rho_bell_adc + K * rho_c * K';
            end
        end
    end
end

disp('Density matrix after ADC:');
disp(simplify(rho_bell_adc));

% Post-noise CNOT operations (same as pre-noise)
rho_final = CNOT_2 * CNOT_1 * rho_bell_adc * CNOT_1' * CNOT_2';

% Measurement operator for escort qubits in state |11⟩
M11 = kron(eye(4), kron([0 0; 0 1], [0 0; 0 1])); 

% Project system onto state where both escort qubits are measured as |1⟩
rho_anc_11 = simplify(M11 * rho_final * M11');
p_11 = trace(rho_anc_11); % Probability of obtaining |11⟩
disp('|11⟩ probability:');
disp(p_11);

% Partial trace over escort qubits to obtain reduced state
reduced_rho = PartialTrace(rho_anc_11, [3, 4]); % Trace out qubits 3 and 4

% Normalize the output state conditioned on successful escort measurement
disp('Output normalized state:');
rout = reduced_rho / p_11;
disp(simplify(rout));
