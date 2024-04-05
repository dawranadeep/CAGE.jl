using Base.Threads, Polynomials, Distances


###################### B-splines #########################
function bspline_recursive1d(x, degree, knots, i)
    # Recursive function to compute B-spline basis functions

    if degree == 0
        B = (x .>= knots[i]) .& (x .< knots[i+1]);
    else
        alpha = (x .- knots[i]) ./ (knots[i+degree] - knots[i]);
        alpha[isinf.(alpha)] .= 0;
        alpha[isnan.(alpha)] .= 0;
        beta = (knots[i+degree+1] .- x) ./ (knots[i+degree+1] - knots[i+1]);
        beta[isinf.(beta)] .= 0;
        beta[isnan.(beta)] .= 0;
        B = alpha .* bspline_recursive1d(x, degree-1, knots, i) .+ beta .* bspline_recursive1d(x, degree-1, knots, i+1);
    end

    return B
end


function bspline_basis_1d(x, degree, knots)
    interrior_knots = knots;
    if interrior_knots[1] == 0
        interrior_knots = interrior_knots[2:end];
    end
    if interrior_knots[end] == 1
        interrior_knots = interrior_knots[1:(end-1)];
    end 
    knots = [zeros(degree+1, 1); interrior_knots; ones(degree+1, 1)]; 
    B = zeros(length(x), length(interrior_knots)+ degree + 1);

    #bspline_recursive1d(x[1:1], degree, knots, 1);

    for i in 1:size(B,2)
        B[:, i] = bspline_recursive1d(x, degree, knots, i);
    end
return B
end




###################### Fourier #########################


function fourier_basis_1d(x, num_args)
    # Compute the Fourier basis functions for a 1D domain

    # Input:
    # x: Coordinates of the evaluation points
    # num_terms: Number of Fourier basis functions to compute

    # Output:
    # B: Matrix of Fourier basis functions evaluated at x

    num_terms = 1 + 2 * num_args;
    B = zeros(length(x), num_terms);

    @threads for i in 1:num_terms
        if i == 1
            B[:, i] = ones(length(x));
        elseif iseven(i)
            B[:, i] = sin.(2π * div(i,2) * x) * sqrt(2);
        else
            B[:, i] = cos.(2π * div(i,2) * x) * sqrt(2);
        end
    end

    return B
end




###################### Legendre #########################

function legendre_basis_1d(x, degree)
    # Compute the Legendre basis functions for a 1D domain

    # Input:
    # x: Coordinates of the evaluation points
    # degree: Degree of Legendre basis functions = number of basis - 2

    # Output:
    # B: Matrix of Legendre basis functions evaluated at x
    P = [];
    if degree==0 
        push!(P, 1);
    else
        push!(P, 1);
        push!(P,[0, 1]);
        for i in 2:(degree+1)            #recurrence relation
            push!(P,  1/i* (  [0; (2*i-1)*P[i]] .- (i-1)* [P[i-1]; 0; 0]));
        end

        for i in 1:(degree+2)           #normalization
            P[i] = P[i] ./ sqrt( 1/(2*i-1) );          # Careful: P[i] = P[i] ./ sqrt( 2/(2*i-1) ) if the range is [-1,1]
        end
    end

    B = zeros(length(x), degree+2);
    @threads for i in 1:(degree+2)
        B[:, i] = Polynomial(P[i]).(2*x .- 1);
    end

    return B
end





###################### Exponential #########################
# exp( - distance/lengthscale)
function exponential_basis_1d(x, lengthscale, knots)
    
    B = exp.( - pairwise(Euclidean(), x, knots) ./ lengthscale);
    B = [ones(length(x), 1) B];
    return B
end




###################### Gaussian #########################
# exp( - distance^2/(2*lengthscale^2))

function gaussian_basis_1d(x, lengthscale, knots)
    
    B = exp.( - pairwise(Euclidean(), x, knots).^2 ./ (2*lengthscale^2) );
    B = [ones(length(x), 1) B];
    return B
end






##################### Collect Everything ###################

function generate_basis_1d(x, basistype = "bspline"; degree=3, knots=[0, 0.5, 1], num_args=10, lengthscale=0.1)
    if basistype == "bspline"
        B = bspline_basis_1d(x, degree, knots);
    elseif basistype == "Fourier"
        B = fourier_basis_1d(x, num_args);
    elseif basistype == "Legendre"
        B = legendre_basis_1d(x, degree);
    elseif basistype == "exponential"
        B = exponential_basis_1d(x, lengthscale, knots);
    elseif basistype == "Gaussian"
        B = gaussian_basis_1d(x, lengthscale, knots);
    else
        throw(DomainError("ErrorMessage: Please read the options of the names of basis function."))
    end
    return B
end
