
###################### Exponential #########################
# exp( - distance/lengthscale)
function exponential_basis_2d(X, lengthscale, knots, intercept)
    
    B = zeros(size(X,1), 0);
    for j in 1:length(lengthscale)
        B = [B  exp.(- 1/ lengthscale[j] * pairwise(Euclidean(), X, knots, dims=1)) ]; 
    end
    if intercept
        B = [ones(size(X,1), 1) B];
    end
    
    return B
end



###################### Gaussian #########################
# exp( - distance^2/2/lengthscale^2)
function Gaussian_basis_2d(X, lengthscale, knots, intercept)

    B = zeros(size(X,1), 0);
    for j in 1:length(lengthscale)
        B = [B  exp.(- 1/(2 * lengthscale[j]^2) * pairwise(Euclidean(), X, knots, dims=1).^2 )];  
    end
    if intercept
        B = [ones(size(X,1), 1) B];
    end

    return B
end




###################### Wendland #########################
function Wendland_basis_2d(X, lengthscale, knots, intercept)
    B = zeros(size(X,1), 0);

    for j in 1:length(lengthscale)
        rr = pairwise(Euclidean(), X, knots, dims=1)/lengthscale[j];
        rr[rr .> 1] .= 1;
        B = [B  (1 .- rr).^6 .* (35* rr.^2 .+ 18 * rr .+ 3)];
    end
    if intercept
        B = [ones(size(X,1), 1) B];
    end

    return B
end





##################### Collect Everything ###################

function generate_basis_2d(X, basistype = "exponential", intercept=false; knots=[0, 0.5, 1], lengthscale=[0.1, 0.3, 0.5])
    if basistype == "exponential"
        B = exponential_basis_2d(X, lengthscale, knots, intercept);
    elseif basistype == "Gaussian"
        B = Gaussian_basis_2d(X, lengthscale, knots, intercept);
    elseif basistype == "Wendland"
        B = Wendland_basis_2d(X, lengthscale, knots, intercept);
    else
        throw(DomainError("ErrorMessage: Please read the options of the names of basis function."))
    end
    return B
end
