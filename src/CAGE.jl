module CAGE;

using Base.Threads, Polynomials, Distances, LinearAlgebra, Trapz, LazyGrids;

include("basis1d.jl");
export bspline_recursive1d, bspline_basis_1d, fourier_basis_1d, legendre_basis_1d, exponential_basis_1d, gaussian_basis_1d, generate_basis_1d;

#include("RowEchelon.jl");
#export rref, rref!;

include("orthonormalize1d.jl");
export orthonormalize_basis_1d;

include("compute_known_KLE_1d.jl");
export compute_known_KLE_1d;


#include("basis2d.jl");
#export bspline_basis_2d;

include("compute_known_MKLE_1d.jl");
export compute_known_MKLE_1d;

include("integrate_simpson_2d.jl");
export integrate_simpson_2d;

end
