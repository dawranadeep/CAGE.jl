include("/Users/dawr2/Desktop/Ranadeep_Daw_Projects/mvcage/github/CAGE/src/CAGE.jl");
using .CAGE, Plots, Distances
using Trapz, LinearAlgebra;
using StatsBase: sample


x = range(0, 1, 100) |> collect;
y = range(0, 1, 200) |> collect;
(xx, yy) = ndgrid(x, y);
X = [xx[:] yy[:]];
lengthscale = 0.1;
knots_x = range(0, 1, 5) |> collect;
knots_y = range(0, 1, 6) |> collect;
(kxx, kyy) = ndgrid(knots_x, knots_y);
knots = [kxx[:] kyy[:]];



B = exponential_basis_2d(X, 0.1, knots );

######### Plot #######
#=
p = plot(x, B[:,1], legend = false);
for i in 2:size(B,2)
    p = plot!(x, B[:,i], legend = false);
end
=#


Q = orthonormalize_basis_2d(X, B);

# Plot randomly 4 of them

idxes = sample(2:size(Q,2), 4; replace=false);
p1 = scatter(x=X[:,1], y=X[:,2], color=Q[:,idxes[1]], legend = false);


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