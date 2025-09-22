
function U = CNOT_on_nqubits(n, control, target)
% Returns the matrix of a CNOT gate acting on an n-qubit system,
% where 'control' and 'target' are indices (1-based) of qubits.

I = eye(2);
X = [0 1; 1 0];

U = zeros(2^n);
for i = 0:(2^n - 1)
    b = dec2bin(i, n) - '0';
    if b(control) == 1
        b(target) = mod(b(target) + 1, 2);
    end
    j = bin2dec(char(b + '0'));
    U(j+1, i+1) = 1;
end
end
