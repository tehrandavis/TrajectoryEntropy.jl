# TrajectoryEntropy

[![Build Status](https://github.com/tehrandavis/TrajectoryEntropy.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/tehrandavis/TrajectoryEntropy.jl/actions/workflows/CI.yml?query=branch%3Amain)

Note: This package is still under development. Please report any issues you encounter. At the moment base functionality is working, but I'm still working on adding some additional features.

This is a Julia package that implements the trajectory entropy analysis method described in Calcagnì, Lombardi, and Sulpizio (2017): "Analyzing spatial data from mouse tracker methodology: An entropic approach" DOI 10.3758/s13428-016-0839-5. 

The code in this repository is a Julia implementation of the MATLAB code provided by the authors of the paper, which can be found [here](http://polorovereto.unitn.it/~antonio.calcagni/emot.html).

## Installation

I'm in the process of registering this package with the Julia package registry. In the meantime, you can install it by copying the contents of this repository to your local machine and running the following in the Julia REPL:

```
julia> include("path/to/directory/TrajectoryEntropy.jl")
```

This should load in all of the necessary functions.

## Usage

The main function of this package is `tea()` (short for "trajectory entropy analysis"). This function takes the follwing arguments:

* `x` - a vector of x-coordinates
* `y` - a vector of y-coordinates
* `tea_options`: A dictionary `Dict()` of options for the analysis. The following options are available:
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
* `plot`: Whether or not to generate plots.
* `verbose`: Whether or not to print out information about the analysis.
* `decompose`: Whether or not to perform entropy decomposition (false = ψ only)




