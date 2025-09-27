"""
AnchoredDensestSubgraph.jl

Main module for the Anchored Densest Subgraph research project.
This module provides access to both ADS (Anchored Densest Subgraph) and 
GADS (Generalized Anchored Densest Subgraph) algorithms, along with 
supporting utilities and competitor algorithms.

# Usage
```julia
using AnchoredDensestSubgraph

# Load a graph
A = readIN("lobster.in", "Example_small/")

# Run ADS algorithm
result = GlobalDensestSubgraph(A)

# Run GADS algorithm  
result = SolveLPDensestSubgraph(A)
```
"""
module AnchoredDensestSubgraph

# Core dependencies
using SparseArrays
using LinearAlgebra
using Base
using Dates
using Printf

# Re-export commonly used functions from Base
export SparseMatrixCSC, sparse, nnz, size, length

# Include and export utility functions
include("utils/Utils.jl")

# Export utility functions
export readIN, readRaw, GetDegree, GetComponentAdjacency, GetOrderByDegreeGraphIndices
export RegisterFunctionStamp, RegisterMemoryItem, Memory_item_GDS, Memory_item_GA, Memory_item_IGA, Memory_item_LA

# Include ADS algorithms (Paper 1)
include("ADS/ADS.jl")

# Export ADS functions
export GlobalDensestSubgraph, LocalAnchoredDensestSubgraph, ImprovedGlobalAnchoredDensestSubgraph
export GlobalAnchoredDensestSubgraph, ProcessGlobalAnchoredDensestSubgraph
export ProcessImprovedGlobalAnchoredDensestSubgraph, ProcessLocalAnchoredDensestSubgraph
export FlowNetAlpha, densestSubgraph

# Include GADS algorithms (Paper 2)  
include("GADS/GADS.jl")

# Export GADS functions
export SolveLPDensestSubgraph, DoSolveLocalADS, SolveLPDensestSubgraphLocal
export SetupLPSolver, DEFAULT_LP_SOLVER, SOLVER_LP_ADSS
export CompareResultSets, ProcessAndOutputAlgorithms

# Include case studies (Universal)
include("case_studies/CaseStudies.jl")

# Export case study functions
export ExportSimpleRs, ImportSimpleRs, ProcessCaseStudy

# Include competitor algorithms
include("competitors/Competitors.jl")

# Include experimental code
include("experiments/Experiments.jl")

# Export competitor functions
export MRW, GreedyL, FlowSeed

# Export experimental functions
export ProcessAlgorithms, DoProcessAlgorithms, GenerateUserInputSet
export BulkProcessAndOutputAlgorithms, DoCollectResults

end # module
