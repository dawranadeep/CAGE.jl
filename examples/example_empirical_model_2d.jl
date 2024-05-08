include("/Users/dawr2/Desktop/Ranadeep_Daw_Projects/mvcage/github/CAGE/src/CAGE.jl");
using .CAGE, PyPlot, Distances, LazyGrids
using Trapz, LinearAlgebra, Distributions, SpecialFunctions;
import Random
Random.seed!(1234)

n = 1000;
X = rand(n,2);
D = pairwise(Euclidean(), X, X, dims=1);

s11 =  sqrt(2);
s22 = 1;
nu11 = 0.4;
nu22 = 0.5;
nu12 = 1/2*(nu11 + nu22);
a11 = 10; #%0.3;
a22 = 15; #%0.5;
a12 = 1.2* max(a11, a22); #0.6; has to be greater than both a11, a22
rhs = a11^nu11 * a22^nu22/ a12^(nu11 + nu22) * gamma(nu12)/gamma(nu12 + 1) / sqrt( gamma(nu11) * gamma(nu22)  ) * sqrt( gamma(nu11 + 1) * gamma(nu22 + 1)  );
rho = rand(Uniform(rhs/2, 3/4 * rhs));

function Matern_cov(D, nu, a)
    tmp =  2^(1 - nu)/ gamma(nu) * (a * D).^nu .* besselk.(nu, a * D);
    tmp[diagind(tmp)] .= 1;
    return tmp
end

C11 = s11^2 * Matern_cov(D, nu11, a11);
C22 = s22^2 * Matern_cov(D, nu22, a22);
C12 = rho * s11 * s22 * Matern_cov(D, nu12, a12);
C = [C11 C12; C12' C22];


nobs = 10000;
Y = rand(MvNormal( zeros(2*n,), C), nobs); 
Y1 = Y[1:n, :];
Y2 = Y[(1+n):end, :];


C11_hat = cov(Y1');
C22_hat = cov(Y2');
C12_hat = cov(Y1', Y2');




## Sample pseudo x and y
px = range(0, 1, 100) |> collect;
py = range(0, 1, 200) |> collect;
(xx, yy) = ndgrid(px, py);
pX = [xx[:] yy[:]];
lengthscale = 0.1;
knots_x = range(0, 1, 20) |> collect;
knots_y = range(0, 1, 30) |> collect;
(kxx, kyy) = ndgrid(knots_x, knots_y);
knots = [kxx[:] kyy[:]];


lengthscale = [0.1, 0.3, 0.5];
B11 = zeros(size(pX,1), 0);

for j in 1:length(lengthscale)
    rr = pairwise(Euclidean(), pX, knots, dims=1)/lengthscale[j];
    rr[rr .> 1] .= 1;
    B11 = [B11  (1 .- rr).^6 .* (35* rr.^2 .+ 18 * rr .+ 3)];
end

Q11 = orthonormalize_basis_2d