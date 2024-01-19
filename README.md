# TrajectoryEntropy


Note: This package is still under development. Please report any issues you encounter. At the moment base functionality is working, but I'm still working on adding some additional features.

This is a Julia package that implements the trajectory entropy analysis method described in Calcagnì, Lombardi, and Sulpizio (2017): "Analyzing spatial data from mouse tracker methodology: An entropic approach" DOI 10.3758/s13428-016-0839-5. 

I've renamed this to Trajectory Entropy Analysis (as opposed to EMoT from the original authors) given you don't need to just do this with mouse trajectories.

The code in this repository is a Julia implementation of the MATLAB code provided by the authors of the paper, which can be found [here](http://polorovereto.unitn.it/~antonio.calcagni/emot.html).

## Installation

I'm in the process of registering this package with the Julia package registry. In the meantime, you can install it by copying the contents of this repository to your local machine and running the following in the Julia REPL:

```
julia> include("path/to/directory/TrajectoryEntropy.jl")
```

This should load in all of the necessary functions.

## Usage

The main function of this package is `tea()`, short for "trajectory entropy analysis". This function takes the following arguments:

* `x` - a vector of x-coordinates
* `y` - a vector of y-coordinates
* `tea_options`: A dictionary `Dict()` of options for the analysis. The following options are available (suggested defaults are shown):
  * "sa_binsize_min" => 12, # minimum bin size for sensitivity analysis
  * "sa_binsize_max" => 20, # maximum bin size for sensitivity analysis
  * "sa_binsize_step" => 2, # step size for sensitivity analysis
  * "verbose" => true, # whether or not to print out information about the analysis
  * "maxInner" => 25000, # maximum number of iterations for inner loop
  * "maxOuter" => 50, # maximum number of iterations for outer loop
  * "maxFunEvals" => 25000, # maximum number of function evaluations
  * "unique_algorithm" => "julia", # algorithm for finding unique points
  * "unique_tol" => 0.01, # tolerance for finding unique points
  * "decompose" => true # whether or not to perform entropy decomposition
* `binHist`: The bin size for the histogram. Can be an integer, "ss" for Shimazaki & Shinomoto, or "sa" for sensitivity analysis.
* `plot`: Whether or not to generate plots of histograms (entropy calculation).
* `verbose`: Whether or not to print out information about the analysis.
* `decompose`: Whether or not to perform entropy decomposition (false = ψ only)

## Note on the options

As you can see above there are a number of crucial options that can be set for the analysis that you may want to adjust provided your data. In particular, the choice of the number of bins (binHist) and the unique tolerance (unique_tol) will have immediate impact on the obtained results. Here, I go into a little more detail about these options.

* `binHist`: The bin size for the histogram. Can be an integer, "ss" for Shimazaki & Shinomoto, or "sa" for sensitivity analysis.
  * If you choose an integer, the histogram will be computed using that number of bins.
  * If you choose "ss", the optimal bin size will be computed using the Shimazaki & Shinomoto method (https://doi.org/10.1162/neco.2007.19.6.1503).
  * If you choose "sa", the optimal bin size will be computed using a sensitivity analysis method. In this case, the "sa_binsize_min", "sa_binsize_max", and "sa_binsize_step" options will be used to set the range of bin sizes to be used in the sensitivity analysis. This method is described in detail in Birgé & Rozenholc (2006): "How many bins should be put in a regular histogram" (https://doi.org/10.1051/ps:2006001).
* `unique_algorithm`: can take one of two values: "julia" or "matlab". A crucial step in the algorithm is reducing the original series of polar coordinates to a series of unique points (i.e., removing duplicate values within a set tolerance). If julia is selected, the `unique()` function from Julia is applied. If "matlab" is selected, then an algortihm analogous to MATLAB's `uniquetol()` is applied. The default is "julia".
* `unique_tol`: The tolerance for finding unique points. The default is 0.01.
* `decompose` (true/false): Whether or not to perform entropy decomposition. If "false" then only the overal entropy value, ψ, is calculated. If set to true, then ψ is further decomposed into ξ (fast movements) and ζ (slow movements) values (see Calcagnì paper). The default is true.
* `plot` (true/false): Whether or not to generate plots. The default is false.
* `verbose` (true/false): Whether or not to print out information about the analysis. The default is true.
* `maxInner`: Entropy decomposition proceeds by a non-linear optimization algorithm (Kullback-Leibler divergence), this option sets the maximum number of iterations that the optimizer is allowed to perform in the inner loop. In the context of the Fminbox optimizer, this would be the maximum number of iterations for the LBFGS algorithm (see Optim.jl: https://julianlsolvers.github.io/Optim.jl/).
* `maxOuter`: sets the maximum number of iterations that the optimizer is allowed to perform in the outer loop. In the context of the Fminbox optimizer, this would be the maximum number of iterations for the box constraint process.
* `maxFunEvals`: sets the maximum number of function evaluations that the optimizer is allowed to perform. This is a stopping criterion for the optimization process. If the optimizer reaches this limit before finding a minimum, it will stop and return the best solution found so far.

## Example
An example of how to perform this analysis is provided in a Jupyter notebook in the examples folder. This folder includes two example trajectories. This notebook can be viewed directly on GitHub, or downloaded and run locally.

## A note on personal practices

Although the sensitivity analysis is implemented in this package, I'm not currently 100% confident in its unsupervised use. My own personal practice at the moment is to manually loop though a range of bin sizes and inspect the resulting decomposition for optimal results. This includes verifying that the obtained decomposition is correct (i.e., that the ξ and ζ values are in the correct range given the obtained ψ). When employing entropy decomposion, I would recommend proceeding with caution.


  
