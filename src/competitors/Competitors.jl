"""
Competitors.jl

Module for competitor algorithm implementations used for comparison with ADS and GADS algorithms.
This module contains both home-coded competitors (CP_*) and external competitors (CX_*).

# Usage
```julia
using AnchoredDensestSubgraph.Competitors

# Use home-coded competitors
result = MRW(graph, seed_nodes)
result = GreedyL(graph, seed_nodes)
result = FlowSeed(graph, seed_nodes)

# Use external competitors
result = ExternalFlowSeed(graph, seed_nodes)
```
"""
module Competitors

# Core dependencies
using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using StatsBase
using Random
using Base

# Include utility functions
include("../utils/Utils.jl")

# Include home-coded competitors
include("home_coded/flowseed.jl")
include("home_coded/greedyl.jl")
include("home_coded/mrw.jl")
include("../experiments/flowseed_compare_test.jl")

# Include external competitors
include("external/flowseed.jl")

# Export competitor functions
export MRW, GreedyL, FlowSeed
export ExternalFlowSeed

# Export comparison and testing functions
export FlowSeedCompareTest

end # module
