include("F:/Now/CAGE.jl/src/CAGE.jl");
using .CAGE, Test

############ Test bspline ###########

n = 5000
X = rand(n,2);
B = 2;