function tea(x, y, tea_options; binHist = "sa", plot = false, verbose = false, decompose = true)
    """
    # TEA: Trajectory Entropy Analysis
    
    ```
    tea(x, y, tea_options; binHist = "sa", plot = false, verbose = false, decompose = true)
    ```

    main function for performing trajectory entropy analysis.

    # Arguments
    - `x`: An array of x-coordinates.
    - `y`: An array of y-coordinates.
    - `tea_options`: A dictionary of options for the analysis.
            Dict(
                "sa_binsize_min" => 12,
                "sa_binsize_max" => 20,
                "sa_binsize_step" => 2,
                "manual_bin_size" => 12,
                "verbose" => true,
                "maxInner" => 25000,
                "maxOuter" => 50,
                "maxFunEvals" => 25000,
                "unique_algorithm" => "julia",
                "unique_tol" => 0.01,
                "decompose" => true
            )
    - `binHist`: The bin size for the histogram. Can be an integer, "ss" for Shimazaki & Shinomoto, or "sa" for sensitivity analysis.
    - `plot`: Whether or not to generate plots.
    - `verbose`: Whether or not to print out information about the analysis.
    - `decompose`: Whether or not to perform entropy decomposition (false = ψ only)

    # Returns
    A dictionary of results from the analysis.
    """
    
    # Set TEA parameters
    value_tolerance = tea_options["unique_tol"]
    maxInner = tea_options["maxInner"]
    maxOuter = tea_options["maxOuter"] 
    maxFunEvals = tea_options["maxFunEvals"]
    sa_Tmin = tea_options["sa_binsize_min"]
    sa_Tmax = tea_options["sa_binsize_max"]


    # Compute the polar coordinates
    polar_degs = coord_to_deg(x,y)    
    rescaled_polar_degs = polar_degs .+ 180
    rescaled_polar_degs[rescaled_polar_degs .> 360] = rescaled_polar_degs[rescaled_polar_degs .> 360] .- 360
    theta = deg2rad.(rescaled_polar_degs)[2:end]    
    
    

    #Find the optimim bin size for the histogram using Shimazaki & Shinomoto
    if binHist == "ss"
        min_bins = sa_Tmin
        max_bins = sa_Tmax
        opt_bin_num = sshist(theta, min_bins, max_bins) 
    elseif binHist == "sa"
        sensitivity_results = sensitivity_binsize(x, y, tea_options)
        opt_bin_num = sensitivity_results["opt_bin_num"]
    else
        opt_bin_num = binHist
    end
    
    if !isa(opt_bin_num, Int)
        opt_bin_num = Int(opt_bin_num)
    end
    

    # setting the bins for the histogram
    xHist = range(0, stop=3.14, length=opt_bin_num+1)

    H_theta = StatsBase.fit(Histogram, theta, xHist).weights
    push!(H_theta, 0)

    # Compute the H_theta0 by getting the unique values of theta
    if tea_options["unique_algorithm"] == "julia"
        theta_0 = round.(theta, digits=3) |> unique
    elseif tea_options["unique_algorithm"] == "matlab"
        theta_0 = uniquetol(theta, value_tolerance)
    else
        println("Invalid unique algorithm")
    end
        
    H_theta0 = StatsBase.fit(Histogram, theta_0, xHist).weights
    push!(H_theta0, 0)


    # entropy decomposition
    entropy_results = entropy_decomposition(H_theta, H_theta0, xHist, maxInner, maxOuter, maxFunEvals; decompose=decompose)
    
    # Plots
    if plot
      # Plot the histogram of theta
      hist_plot = histogram(theta, bins = opt_bin_num, xlabel = "theta", ylabel = "Frequency", title = "θ Values", legend=false)
    else
        hist_plot = 0
    end
    
    # Create a dictionary combining important results / parameters with entropy_decomposition results
    output_results = Dict{String, Any}("hist_plot" => hist_plot, "theta" => theta, "H_theta" => H_theta, "H_theta0" => H_theta0, "opt_bin_num" => opt_bin_num)
    
    for (key, value) in entropy_results
        output_results[key] = value
    end
    
    return output_results
end

