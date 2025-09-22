
function P = projector_on_escort_111(escort_indices, n)
% Returns the projection operator on escort qubits = |111⟩ in an n-qubit system.
% escort_indices: array with indices (1-based) of the escort qubits

P = 1;
for i = 1:n
    if ismember(i, escort_indices)
        P = kron(P, [0 0; 0 1]);  % projector on |1⟩
    else
        P = kron(P, eye(2));
    end
end
end
