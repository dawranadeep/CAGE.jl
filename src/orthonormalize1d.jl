

function orthonormalize_basis_1d(x, B)   

    Q = zeros(size(B));
    
    # Loop through columns
    for i in 1:size(B,2)
        v = B[:, i];
        for j in 1:i-1
            # Projection of current vector onto previous ones
            proj =  trapz(x, v.* Q[:, j])/trapz(x, Q[:,j].* Q[:, j]) * Q[:, j];
            # Subtract projection
            v -= proj
        end;

        nrm = sqrt(trapz(x, v.*v));
        if nrm > eps() 
            Q[:, i] = v/ nrm;
        end
    end
    return Q

end






