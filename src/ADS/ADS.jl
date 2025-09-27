"""
ADS.jl

Module for Anchored Densest Subgraph (ADS) algorithms from Paper 1.
This module contains the flow network-based algorithms for finding
anchored densest subgraphs in graphs.

# Usage
```julia
using AnchoredDensestSubgraph.ADS

# Load a graph
A = readIN("lobster.in", "Example_small/")

# Run global densest subgraph algorithm
result = GlobalDensestSubgraph(A)

# Run local anchored densest subgraph algorithm
result = LocalAnchoredDensestSubgraph(A, [1, 2])
```
"""
module ADS

# Core dependencies
using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using Base

# Include utility functions
include("../Utils_io.jl")
include("../Utils_graph.jl")
include("../Utils.jl")
include("../Memory_tracker.jl")

# Include ADS algorithms
include("algorithms.jl")
include("flow_network.jl")

# Include experiments
include("experiments/query_tests.jl")
include("experiments/degeneracy_tests.jl")

# Export ADS functions
export GlobalDensestSubgraph, LocalAnchoredDensestSubgraph, ImprovedGlobalAnchoredDensestSubgraph
export GlobalAnchoredDensestSubgraph, ProcessGlobalAnchoredDensestSubgraph
export ProcessImprovedGlobalAnchoredDensestSubgraph, ProcessLocalAnchoredDensestSubgraph
export FlowNetAlpha, densestSubgraph

# Export experimental functions
export GenerateUserInputSet, RandomGenerateUserInputSet, DoProcessAlgorithms
export GetSoleAlgorithmIndex, DoOutputPerformanceReports

# Export memory tracking constants
export Memory_item_GDS, Memory_item_GA, Memory_item_IGA, Memory_item_LA

end # module
