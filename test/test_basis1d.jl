include("/Users/dawr2/Desktop/Ranadeep_Daw_Projects/mvcage/github/CAGE/src/CAGE.jl");
using .CAGE, Test

############ Test bspline ###########

x = range(0, 1, 1000) |> collect;
knots = range(0, 1, 10) |> collect;
degree = 3;
B = bspline_basis_1d(x, degree, knots);

@testset "basis Test" begin
    @test [minimum(x), maximum(x)] == [0, 1] || "Domain not within 0 and 1"
    @test (minimum(knots) >= 0) & (maximum(knots) <= 1) || "Domain not within 0 and 1"
    @test sum(isnan.(B)) == 0 || "NA values in bspline"
    #@test size(B) == (length(x), length(knots) + degree - 1) || "Check the dimension of bspline"
end
