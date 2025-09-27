"""
Experiments.jl

Module for experimental code and testing functionality.
This module contains experimental algorithms, comparison tests, and evaluation code.

# Usage
```julia
using AnchoredDensestSubgraph.Experiments

# ADS experiments
result = run_ads_query_test(graph, query_nodes)

# GADS experiments  
result = run_gads_lp_compare_test(graph, parameters)

# FlowSeed comparison
result = run_flowseed_compare_test(graph, seeds)
```
"""
module Experiments

# Core dependencies
using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using StatsBase
using Random
using Base

# Include experimental modules
include("ads_query_tests.jl")
include("ads_degeneracy_tests.jl")
include("gads_lp_compare_test.jl")
include("flowseed_compare_test.jl")

# Export experimental functions
export run_ads_query_test, run_ads_degeneracy_test
export run_gads_lp_compare_test
export run_flowseed_compare_test

end # module


