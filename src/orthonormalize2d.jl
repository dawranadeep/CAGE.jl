

function orthonormalize_basis_2d(X, B)   
    B2 = reshape(B, (size(x,1), size(y,1), size(B,2)));
    Q2 = zeros(size(B2));
    x1 = unique(X[:,1]) |> sort;
    x2 = unique(X[:,2]) |> sort;


    
    # Loop through columns
    for i in 1:size(B2,3)
        v = B2[:, :, i];
        for j in 1:i-1
            # Projection of current vector onto previous ones
            proj =  trapz(x1, trapz(x2, v.* Q2[:, :, j]))/trapz(x1, trapz(x2, Q2[:, :, j].* Q2[:, :, j])) * Q2[:, :, j];
            # Subtract projection
            v -= proj
        end;

        nrm = sqrt(trapz(x1, trapz(x2, v.*v)));
        if nrm > eps() 
            Q2[:, :, i] = v/ nrm;
        end
    end
    return reshape(Q2, (size(x,1) * size(y,1), size(B,2)));

end






