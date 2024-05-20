module rdtest

__precompile__()

using LinearAlgebra, LazyGrids, Trapz, Distances;

include("simpmat2d.jl");
export simpmat2d;

include("basis2d.jl");
export exponential_basis_2d, Gaussian_basis_2d, Wendland_basis_2d, generate_basis_2d;

include("orthogonalize_basis2d.jl");
export orthogonalize_basis2d_v1, orthogonalize_basis2d_v2;

end
