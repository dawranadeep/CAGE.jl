
function compute_known_KLE_1d(x, C, B)
    N_terms = size(B, 2);
    A = zeros(N_terms, N_terms);
    for i in 1:N_terms
        for j in i:N_terms
            A[i, j] = first(trapz(x, mapslices(t -> trapz(x, t), C .* (B[:,i] * B[:,j]'), dims=1)));
            A[j, i] = A[i, j];
        end
    end
    
    EE = eigen(A);
    Psi = B * EE.vectors;
    lam = EE.values;

    return Psi, lam

end

