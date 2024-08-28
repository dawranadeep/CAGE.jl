include("F:/Now/CAGE.jl/src/CAGE.jl");
using .CAGE, Test
using LazyGrids, Statistics


############ Test bspline ###########

n = 5000
X = rand(n,2);

knots_x = range(0, 1, 20) |> collect;
knots_y = range(0, 1, 30) |> collect;
(kxx, kyy) = ndgrid(knots_x, knots_y);
knots = [kxx[:] kyy[:]];

B = exponential_basis_2d(X, 0.1, knots, false );


