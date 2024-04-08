include("/Users/dawr2/Desktop/Ranadeep_Daw_Projects/mvcage/github/CAGE/src/CAGE.jl");
using .CAGE, Plots, Distances
using Trapz, LinearAlgebra;

x = range(0, 1, 1000) |> collect;
knots = range(0, 1, 25) |> collect;
degree = 3; # For bspline
num_args = 3; # For Fourier
lengthscale = 0.1;

B = generate_basis_1d(x, "Gaussian"; knots=knots );
#B = generate_basis_1d(x, "Legendre"; degree = 20);

######### Plot #######
#=
p = plot(x, B[:,1], legend = false);
for i in 2:size(B,2)
    p = plot!(x, B[:,i], legend = false);
end
=#


Q = orthonormalize_basis_1d(x, B);
p = plot(x, Q[:,1], legend = false);
for i in 2:size(Q,2)
    p = plot!(x, Q[:,i], legend = false);
end


W2 = zeros(size(Q,2), size(Q,2));
for j in 1:size(Q,2)
    for k in j:size(Q,2)
        W2[j,k] = trapz(x, Q[:, j] .* Q[:, k]);
        W2[k,j] = W2[j, k];
    end
end


round.(W2, digits=4)

########        Test KLE ##########


C = exp.( -pairwise(Euclidean(), x, x)/5);
(psi, ee) = compute_known_KLE_1d(x, C, Q);


W2 = zeros(size(Q,2), size(Q,2));
for j in 1:size(Q,2)
    for k in j:size(Q,2)
        W2[j,k] = trapz(x, psi[:, j] .* psi[:, k]);
        W2[k,j] = W2[j, k];
    end
end

p = plot(x, psi[:,1], legend = false);
for i in 2:size(psi,2)
    p = plot!(x, psi[:,i], legend = false);
end