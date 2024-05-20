include("/Users/dawr2/Desktop/Ranadeep_Daw_Projects/mvcage/github/CAGE/src/CAGE.jl");
using .CAGE
using LazyGrids
# using StatsBase: sample


x = range(0, 1, 101) |> collect;
y = range(0, 1, 201) |> collect;
(xx, yy) = ndgrid(x, y);
X = [xx[:] yy[:]];

