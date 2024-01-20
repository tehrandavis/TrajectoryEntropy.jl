function sshist(x, min_bins, max_bins)
    """
    sshist(x, min_bins, max_bins)
    

    Function `sshist' returns the optimal number of bins in a histogram
    used for density estimation.
    
    Optimization principle is to minimize expected L2 loss function between 
    the histogram and an unknown underlying density function.
    An assumption made is merely that samples are drawn from the density
    independently each other.

    The optimal binwidth D* is obtained as a minimizer of the formula, 
    (2K-V) / D^2,
    where K and V are mean and variance of sample counts across bins with width D.
    Optimal number of bins is given as (max(x) - min(x)) / D*.

    For more information, visit 
    http://2000.jukuin.keio.ac.jp/shimazaki/res/histogram.html

    Original paper:
    Hideaki Shimazaki and Shigeru Shinomoto
    A method for selecting the bin size of a time histogram
    Neural Computation 19(6), 1503-1527, 2007
    http://dx.doi.org/10.1162/neco.2007.19.6.1503

    # Arguments
    - `x`: An array of values.
    - `min_bins`: The minimum number of bins.
    - `max_bins`: The maximum number of bins.
    
    # Returns
    The optimal number of bins.
    
    """
    x_min, x_max = minimum(x), maximum(x)

    N_MIN = min_bins  # Min number of bins (integer), must be > 1
    N_MAX = max_bins # Max number of bins (integer)
    SN = 30 # of partitioning positions for shift average
        
    N = N_MIN:N_MAX # number of bins
    D = (x_max-x_min)./N        # Bin size vector
    C = zeros(length(D), SN)

    # Computation of the cost function
    for i in eachindex(N)
        shift = range(0, stop=D[i], length=SN)
        
        for p in 1:SN
            edges = range(x_min + shift[p] - D[i]/2, stop=x_max + shift[p] - D[i]/2, length=N[i]+1) # Bin edges
            ki = fit(Histogram, x, edges).weights     # Count number of events in bins
            k = mean(ki)                         # Mean of event count
            v = sum((ki .- k).^2)/N[i]              # Variance of event count
            C[i, p] = (2*k-v)/(D[i]^2)              # Cost Function
        end
    end

    # Optimal bin size Selection
    cmin = minimum(C)
    idx  = findfirst(C .== cmin)
    optN = N[idx[1]]
    
    if optN % 2 != 0
        optN = optN - 1
    end
    
    return optN
end



##---- Sensitivity Analysis ----##

function sensitivity_binsize(x,y, tea_options)
    # Set TEA parameters for the sensitivity analysis
    bin_min = tea_options["sa_binsize_min"]
    bin_max = tea_options["sa_binsize_max"]
    bin_step = tea_options["sa_binsize_step"]
    verbose = tea_options["verbose"]

    # Run TEA over different T
    if verbose println("@@ Optimizing number of bins") end
    c = 0
    countNaN = 0
    D = []
    for bin_num = bin_min:bin_step:bin_max
        if verbose println(bin_num) end
        c += 1
        #binHist = k
        try
            res = tea(x, y, tea_options; binHist = bin_num, plot=false, verbose=true)
            if all(x -> haskey(res, x), ["psi", "csi", "zeta1", "zeta2"])
                D = vcat(D, [bin_num 0 res["psi"] res["csi"] res["zeta1"] res["zeta2"]])
            else
                println("Missing keys in the result: ", setdiff(["psi", "csi", "zeta1", "zeta2"], keys(res)))
            end
        catch e
            println("Error in iteration $bin_num: $e")
        end
        
        # to do, build in exception handling for when the optimization fails
    end

    # Compute best bin size
    C = [(round.(D[:,5], digits=2) .> round.(D[:,6], digits=2)) (round.(D[:,5], digits=2) .< round.(D[:,6], digits=2)) (round.(D[:,5], digits=2) .== round.(D[:,6], digits=2))]
    _, f = findmax(sum(C, dims=1))
    iid = C[:, f[2]]
    Dd = D[iid, :]  

    tBest = zeros(3)
    if median(Dd[:,4]) > 0.1
        tBest[1] = maximum(D[findall(abs.(Dd[:,4] .- median(Dd[:,4])) .== minimum(abs.(Dd[:,4] .- median(Dd[:,4])))), 1])
    end
    if median(Dd[:,5]) > 0.1
        tBest[2] = maximum(D[findall(abs.(Dd[:,5] .- median(Dd[:,5])) .== minimum(abs.(Dd[:,5] .- median(Dd[:,5])))), 1])
    end
    if median(Dd[:,6]) > 0.1
        tBest[3] = maximum(D[findall(abs.(Dd[:,6] .- median(Dd[:,6])) .== minimum(abs.(Dd[:,6] .- median(Dd[:,6])))), 1])
    end

    tBest = tBest[tBest .> 0]

    if !isempty(tBest)
        tOpt = round(mean(tBest))
        if tOpt % 2 != 0
            tOpt = tOpt - 1
        end
        accOpt = (sum(C[:,f[2]]) / length(C[:,f[2]])) * 100  # compute accuracy of the best binsize
    else
        try
            tOpt = maximum(D[findmin(abs.(Dd[:,3] .- median(Dd[:,3]))), 1])  # with respect to: PSI
        catch
            accOpt = NaN
        end
    end

    if verbose
        println("done: $tOpt")
        println("($accOpt%)")
        println()
    end

    # Save results
    results = Dict("data" => D, "opt_bin_num" => tOpt, "acc" => accOpt)
end