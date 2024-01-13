function uniquetol(A, tol)
    
    """
    uniquetol(A, tol)
    
    Returns the unique elements in A using tolerance tol. 
    Two values, u and v, are within tolerance if:
    abs(u-v) <= tol*max(abs(A(:))) 
    
    That is, uniquetol scales the tol input based on the magnitude of the data.
    
    # Arguments
    - `A`: An array of values.
    - `tol`: A tolerance value.
    
    # Returns
    An array of unique values.
    """
    
    tol_scaled = tol * maximum(abs.(A))  # scale tolerance
    sort!(A)  # sort A in place
    C = [A[1]]  # initialize output with the first element of A

    for i in eachindex(A)[2:end]
        if !(A[i] - C[end] <= tol_scaled)
            push!(C, A[i])
        end
    end

    return C
end