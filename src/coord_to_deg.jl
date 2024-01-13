function coord_to_deg(x,y)
    """
    coord_to_deg(x, y)

    Compute the angle in degrees between the x-axis and the vector (x, y).

    # Arguments
    - `x`: An array of x-coordinates.
    - `y`: An array of y-coordinates.

    # Returns
    An array of angles in degrees.
    """
    p = atan.(y, x) .* -180 / π
    p .= ifelse.(p .< 0, p .+ 360, p)
    return p
    
        # p = similar(x)
        # for i in eachindex(x)
        #     deg = -((atan(y[i], x[i]))*360)/(π*2)
        #     if deg > 0 && deg < 180
        #         p[i] = deg
        #     elseif deg == -180 
        #         p[i] = deg*-1
        #     elseif deg == 0.0 
        #         p[i] = 360.00
        #     elseif deg < 0.0 && deg > -180 
        #         p[i] = 360.00-(deg*-1)
        #     end
        # end
        # return p
end