module CAGE;

include("basis1d.jl");
export bspline_recursive1d, bspline_basis_1d, fourier_basis_1d, legendre_basis_1d, exponential_basis_1d, gaussian_basis_1d, generate_basis_1d;

#include("RowEchelon.jl");
#export rref, rref!;

include("orthonormalize1d.jl");
export orthonormalize_basis_1d;

#include("orthonormalize.jl");
#export orthonormalized_basis_1d;
end
