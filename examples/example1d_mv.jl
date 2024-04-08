include("/Users/dawr2/Desktop/Ranadeep_Daw_Projects/mvcage/github/CAGE/src/CAGE.jl");
using .CAGE, Plots, Distances, LazyGrids
using Trapz, LinearAlgebra;








########        Test MKLE ##########

x = range(0, 1, 1000) |> collect;
C11 = exp.( -pairwise(Euclidean(), x, x)/5);
C22 = exp.( -pairwise(Euclidean(), x, x)/5);
C12 = 0.0 * C11;

B1 = generate_basis_1d(x, "Fourier"; num_args = 10);
B2 = generate_basis_1d(x, "Fourier"; num_args = 10);

Q1 = orthonormalize_basis_1d(x, B1);
Q2 = orthonormalize_basis_1d(x, B2);


(psi12, lam12, psi11, lam11, psi22, lam22) = compute_known_MKLE_1d(x, C11, C22, C12, Q1, Q2);

###### Check eigenfunctions ########
tmp = [psi12[1] psi12[2]];
W2 = zeros(size(tmp,2), size(tmp,2));
for i in 1:size(tmp, 2)
    for j in 1:size(tmp, 2)
        W2[i,j] = trapz(x, tmp[:,i] .* tmp[:, j] );
        W2[j,i] = W2[i,j];
    end
end



p = plot(x, tmp[:,1], legend = false);
for i in 2:size(tmp,2)
    p = plot!(x, tmp[:,i], legend = false);
end