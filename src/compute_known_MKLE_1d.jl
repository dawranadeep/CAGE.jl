
function compute_known_MKLE_1d(x, C11, C22, C12, B1, B2)
    
    (tmpPsi1, tmplam1) = compute_known_KLE_1d(x, C11, B1);
    (tmpPsi2, tmplam2) = compute_known_KLE_1d(x, C22, B2);


    N_terms1 = size(B1, 2);
    N_terms2 = size(B2, 2);
    
    A12 = zeros(N_terms1, N_terms2);
    for i in 1:N_terms1
        for j in 1:N_terms2
            #print(i, "\t", j, "\n");
            A12[i, j] = first(trapz(x, mapslices(t -> trapz(x, t), C12 .* (tmpPsi1[:,i] * tmpPsi2[:,j]'), dims=1)));
            #A12[j, i] = A12[i, j];
        end
    end
    
    AA = [diagm(tmplam1) A12; A12' diagm(tmplam2)];
    EE = eigen(AA);

    # Remove negative eigenvalues
    ii = EE.values .> 0;
    
    Psi12 = [];
    push!(Psi12, B1 * EE.vectors[1:N_terms1, ii]);
    push!(Psi12, B2 * EE.vectors[(1+N_terms1):end, ii]);
    
    lam12 = EE.values[ii];

    return Psi12, lam12, tmpPsi1, tmplam1, tmpPsi2, tmplam2
end

