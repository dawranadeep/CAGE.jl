function orthogonalize_basis2d_simpson(x1_l, x1_u, x2_l, x2_u, nx1, nx2, B; knots=[0, 0.5, 1], lengthscale=[0.1, 0.3, 0.5])
    
    # Compute the Simpson's matrix for integration
    XW = simpmat2d(x1_l, x1_u, x2_l, x2_u, nx1, nx2);
    sW = XW[:, end];
    
    # Compute the basis functions
    pB  = generate_basis_2d(XW[:, 1:2], "Gaussian", false, lengthscale=[0.5, 0.3, 0.1], knots=knots );
    

    

    # Compute orthogonal matrix and Gram-Schmidt post multiplication
    W = [];
    qB = zeros(size(pB));

    for i in 1:size(pB,2)
        v = pB[:, i];
        for j in 1:(i-1)
            # Projection of current vector onto previous ones
            proj = sum(v .* qB[:, j] .* sW) ;
            # Subtract projection
            v = v - proj * qB[:, j];
        end;

        nrm = sqrt( sum(v .* v .* sW) );
        if nrm > eps() 
            qB[:, i] = v/ nrm;
        end;
        
    end




    
    
end








function orthogonalize_basis2d_v1(x1_l, x1_u, x2_l, x2_u, nx1, nx2, basistype, intercept, threshold; knots=[0, 0.5, 1], lengthscale=[0.1, 0.3, 0.5])
    

    grid_x1 = range(x1_l, x1_u, nx1+1);
    grid_x2 = range(x2_l, x2_u, nx2+1);
    xx1, xx2 = ndgrid(grid_x1, grid_x2);
    XX = [xx1[:] xx2[:]];
    
    # Compute the basis functions
    pB  = generate_basis_2d(XX, basistype, intercept, lengthscale=lengthscale, knots=knots );
    pB  = reshape(pB, (nx1+1, nx2+1, size(pB,2)) );

    

    # Compute orthogonal matrix and Gram-Schmidt post multiplication
    W_proj = []; W_nrm = [];
    qB = zeros(size(pB));

    # Loop through columns
    for i in 1:size(pB,3)
        
        v = pB[:, :, i];
        
        this_proj = [];
        for j in 1:(i-1)
            # Projection of current vector onto previous ones
            proj =  trapz(grid_x1, trapz(grid_x2, pB[:, :, i].* qB[:, :, j]));
            # Subtract projection
            v = v - proj * qB[:, :, j];
            push!(this_proj, proj);
        end;
        push!(W_proj, this_proj);

        nrm = sqrt(trapz(grid_x1, trapz(grid_x2, v.*v)));
        if nrm < threshold
            nrm = Inf;
        end;
        qB[:, :, i] = v/ nrm;
        push!(W_nrm, nrm);
        
    end

    qB = reshape(qB, ( (nx1+1) * (nx2+1), size(qB,3)));

    return qB, W_proj, W_nrm
    
    
end








function orthogonalize_basis2d_v2(x1_l, x1_u, x2_l, x2_u, nx1, nx2, basistype, intercept, threshold; knots=[0, 0.5, 1], lengthscale=[0.1, 0.3, 0.5])
####### Version 2
####### First orthogonalize, then normalize    

    grid_x1 = range(x1_l, x1_u, nx1+1);
    grid_x2 = range(x2_l, x2_u, nx2+1);
    xx1, xx2 = ndgrid(grid_x1, grid_x2);
    XX = [xx1[:] xx2[:]];
    
    # Compute the basis functions
    pB  = generate_basis_2d(XX, basistype, intercept, lengthscale=lengthscale, knots=knots );
    pB  = reshape(pB, (nx1+1, nx2+1, size(pB,2)) );

    

    # Compute orthogonal matrix and Gram-Schmidt post multiplication
    W_proj = []; W_nrm = [];
    qB = zeros(size(pB));

    # Loop through columns
    for i in 1:size(pB,3)
        
        v = pB[:, :, i];
        
        this_proj = [];
        for j in 1:(i-1)
            # Projection of current vector onto previous ones
            proj =  trapz(grid_x1, trapz(grid_x2, pB[:, :, i].* qB[:, :, j]))/trapz(grid_x1, trapz(grid_x2, qB[:, :, j].^2));
            # Subtract projection
            v = v - proj * qB[:, :, j];
            push!(this_proj, proj);
        end;
        push!(W_proj, this_proj);

        qB[:, :, i] = v;
        nrm = sqrt(trapz(grid_x1, trapz(grid_x2, v.*v)));
        if nrm < threshold
            nrm = Inf;
        end; 
        qB[:, :, i] = v/ nrm;
        push!(W_nrm, nrm);
        
    end

    
    

    qB = reshape(qB, ( (nx1+1) * (nx2+1), size(qB,3)));

    return qB, W_proj, W_nrm
    
    
end
