

function bspline_recursive1d(x, degree, knots, i)
    # Recursive function to compute B-spline basis functions

    if degree == 0
        B = (x .>= knots[i]) .& (x .< knots[i+1])
    else
        alpha = (x .- knots[i]) ./ (knots[i+degree] - knots[i]);
        alpha[isinf.(alpha)] .= 0;
        beta = (knots[i+degree+1] .- x) ./ (knots[i+degree+1] - knots[i+1]);
        beta[isinf.(beta)] .= 0;
        B = alpha .* bspline_recursive1d(x, degree-1, knots, i) .+ beta .* bspline_recursive1d(x, degree-1, knots, i+1);
    end

    return B
end


function bspline_basis_1d(x, degree, knots)
    tmp_n = length(knots) - (degree + 1);
    B = zeros(length(x), tmp_n);

    Threads.@threads for i in 1:tmp_n
        B[:, i] = bspline_recursive(x, degree, knots, i)
    end
return B
end


