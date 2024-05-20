function simpmat2d(x1_l, x1_u, x2_l, x2_u, nx1, nx2)
    # Ranadeep Daw
    # Courtesy: Whayne Padden's matlabexchange code

    # Only calculates the Simpson's matrix for 2d. Does not integrate.
    # x1_l (lower) and x1_u (upper) are limits of first axis of integration
    # x2_l (lower) and x2_u (upper) are limits of second axis of integration
    # nx1 and nx2 are number of intervals in the first and second axis, respectively.
    # nx1, nx2 should be even. If not, the program forces it.

    # Example: To integrate a function f(x,y) over the region 0 <= x <= 4, 8 <= y <= 12,
    # use simpmat2d(0, 4, 8, 12, 100, 200)
    # where 100 and 200 intervals are used in the x and y direction 


    # Ensure the numbers of intervals are even
    nx1 = Int(2*ceil(nx1/2));
    nx2 = Int(2*ceil(nx2/2));

    # Calculate odd and even indices. Note: usual notation starts from 0,
    # so odd numbers are 2:2:nx1
    io1 = 2:2:nx1;
    ie1 = 3:2:nx1-1;
    io2 = 2:2:nx2;
    ie2 = 3:2:nx2-1;

    # Contributions
    # End points: 1; Odds : 4; Evens : 2 
    w1 = ones(nx1+1);
    w1[io1] .= 4;
    w1[ie1] .= 2;
    w2 = ones(nx2+1);
    w2[io2] .= 4;
    w2[ie2] .= 2;
    w = w1 .* w2';

    hx = (x1_u - x1_l)/nx1;
    hy = (x2_u - x2_l)/nx2;
    SW = w * hx * hy/ 9.0;

    grid_x1 = range(x1_l, x1_u, nx1+1);
    grid_x2 = range(x2_l, x2_u, nx2+1);
    xx1, xx2 = ndgrid(grid_x1, grid_x2);
    XW = [xx1[:] xx2[:] SW[:]];
    
    return(XW)
end
