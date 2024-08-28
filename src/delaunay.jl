
function delaunay_integration(X, Y)
# Calculate a discrete integration with locations (2d) X and observations (1d) Y
# Step 1: Triangulate using Delaunay triangulation on X
# Step 2: Multiply areas with average observations
    mesh = delaunay(X);


    function area(a, b, c)
        return abs( 0.5 * ( a[1] * (b[2] − c[2]) + b[1] * (c[2] − a[2]) + c[1] * (a[2] − b[2]) ))
    end

    function rowwise_area(X, mesh, t)
        return area(X[mesh.simplices[t,1], :], X[mesh.simplices[t,2], :], X[mesh.simplices[t,3], :])
    end


    areas = zeros(size(mesh.simplices, 1), 1);
    for j in 1:size(mesh.simplices, 1)
        areas[j] =     rowwise_area(X, mesh, j);
    end


    
    avgYs = mapslices(t -> mean(t), Y[mesh.simplices], dims=2);
    integral = areas .* avgYs |> sum;

    return integral

end