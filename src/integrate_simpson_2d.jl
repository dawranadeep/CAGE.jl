
function integrate_simpson_2d(x1g, x2g, gridy)
    # Locations should be a meshgrid of size n1 X n2.
    # gridx1 and gridx2 should be NX X NY grid of the first and second location parameters
    # gridy should be the gridded observvvation

    
    xxg, yyg = ndgrid(x1g,x2g);

    return xxg

end

