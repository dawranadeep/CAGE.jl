

function orthonormalize_basis_2d(X, B; B2)   
    B_tmp = reshape(B, (size(x,1), size(y,1), size(B,2)));
    Q_tmp = zeros(size(B_tmp));
    x1 = unique(X[:,1]) |> sort;
    x2 = unique(X[:,2]) |> sort;

    
    
    # Loop through columns
    for i in 1:size(B_tmp,3)
        v = B_tmp[:, :, i];
        for j in 1:i-1
            # Projection of current vector onto previous ones
            proj =  trapz(x1, trapz(x2, v.* Q_tmp[:, :, j]))/trapz(x1, trapz(x2, Q_tmp[:, :, j].* Q_tmp[:, :, j])) * Q_tmp[:, :, j];
            # Subtract projection
            v -= proj
        end;

        nrm = sqrt(trapz(x1, trapz(x2, v.*v)));
        if nrm > eps() 
            Q_tmp[:, :, i] = v/ nrm;
        end
    end
    return reshape(Q_tmp, (size(x,1) * size(y,1), size(B,2)));

end






