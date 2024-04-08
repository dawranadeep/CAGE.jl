
function bspline_basis_2d(x, y, degree, knotsx, knotsy)
    @assert ndims(x) == 1 & ndims(y) == 1;
    @assert ndims(degree) <= 1;
    @assert ndims(knotsx) == 1 & ndims(knotsy) == 1;
    if ndims(degree) == 0
        degree = [degree, degree];
    end;
    xB = bspline_basis_1d(x, degree[1], knotsx);
    yB = bspline_basis_1d(y, degree[2], knotsy);
    xyB = zeros(size(x,1), size(y,1), size(xB,2)* size(yB,2));
    for j in 1:size(xyB,3)
        #print(j, "\n")
        j1 = Int(1 + floor((j - 1)/size(yB,2)));
        j2 = 1 + (j - 1) % size(yB,2);
        xyB[:, :, j] = xB[:, j1] * yB[:, j2]';
    end

    return xyB
end





###################### Exponential #########################
# exp( - distance/lengthscale)
function exponential_basis_1d(X, lengthscale, knots)
    
    if ndims(lengthscale) == 0
        B = exp.( - pairwise(Euclidean(), X, knots) ./ lengthscale);
        B = [ones(length(x), 1) B];
    end
    return B
end




###################### Gaussian #########################
# exp( - distance^2/(2*lengthscale^2))

function gaussian_basis_1d(x, lengthscale, knots)
    
    B = exp.( - pairwise(Euclidean(), x, knots).^2 ./ (2*lengthscale^2) );
    B = [ones(length(x), 1) B];
    return B
end




