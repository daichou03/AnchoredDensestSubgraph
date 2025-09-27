"""
Utils.jl

Module for shared utility functions used across ADS, GADS, case studies, and competitors.
This module contains core utilities, graph operations, I/O functions, memory tracking,
and result collection functionality.

# Usage
```julia
using AnchoredDensestSubgraph.Utils

# Core utilities
result = some_core_function()

# Graph operations
graph_info = analyze_graph(adjacency_matrix)

# I/O operations
data = read_graph_file("graph.in")

# Memory tracking
track_memory_usage("operation_name")

# Result collection
collect_experiment_results(results)
```
"""
module Utils

# Core dependencies
using SparseArrays
using MAT
using LinearAlgebra
using StatsBase
using Random
using Base
using CSV
using DataFrames

# Include utility modules
include("core.jl")
include("graph.jl")
include("io.jl")
include("warmup.jl")
include("memory_tracker.jl")
include("collect_results.jl")

# Export core utility functions
export readIN, readRaw, GetDegree, GetComponentAdjacency, GetOrderByDegreeGraphIndices
export RegisterFunctionStamp, RegisterMemoryItem, Memory_item_GDS, Memory_item_GA, Memory_item_IGA, Memory_item_LA

# Export graph utility functions
export toTransitionGraph, GetAdjacency, GetVolume, GetConductance
export GetRefinedSet, GetRefinedSetFromR, GenerateReferenceSetFixedWalks

# Export I/O utility functions
export folderString, emptyStringArray, readAnchors, readCompsets
export writeResults, writePerformanceReport

# Export memory tracking functions
export RegisterFunctionStamp, RegisterMemoryItem
export Memory_item_GDS, Memory_item_GA, Memory_item_IGA, Memory_item_LA

# Export result collection functions
export ProcessAlgorithms, DoProcessAlgorithms, GenerateUserInputSet
export BulkProcessAndOutputAlgorithms, DoCollectResults

# Export warmup functions
export WarmupCoreAlgorithms

end # module

