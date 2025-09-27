"""
GADS.jl

Module for Generalized Anchored Densest Subgraph (GADS) algorithms from Paper 2.
This module contains the linear programming-based algorithms for finding
generalized anchored densest subgraphs in graphs.

# Usage
```julia
using AnchoredDensestSubgraph.GADS

# Load a graph
A = readIN("lobster.in", "Example_small/")

# Run LP densest subgraph algorithm
result = SolveLPDensestSubgraph(A)

# Run local anchored densest subgraph with weights
R = [1, 2]
weight = [2, 1, 0, 0, 0, 0, -1]
result = DoSolveLocalADS(SOLVER_LP_ADSS, A, R, false, false, DEFAULT_LP_SOLVER, weight)
```
"""
module GADS

# Core dependencies
using SparseArrays
using MAT
using LinearAlgebra
using Base
using JuMP
using CSV
using DataFrames
using StatsBase

# Include utility functions
include("../utils/Utils.jl")

# Include GADS algorithms
include("lp_consts.jl")
include("lp_load_solver.jl")
include("algorithms.jl")
include("lp_evaluation.jl")

# Include experiments
include("../experiments/gads_lp_compare_test.jl")

# Export GADS functions
export SolveLPDensestSubgraph, DoSolveLocalADS, SolveLPDensestSubgraphLocal
export SetupLPSolver, DEFAULT_LP_SOLVER, SOLVER_LP_ADSS
export CompareResultSets, ProcessAndOutputAlgorithms
export BulkProcessAndOutputAlgorithms, ProcessAndOutputLPFixedSizes

# Export constants
export LP_ALGORITHM_NAMES, LP_COMP_RESULT_NAMES, LP_EVAL_RESULT_NAMES
export STATS_NAMES, EVAL_NAMES, FOLDER_LP_COMP_RESULTS, FOLDER_LP_EVAL_RESULTS

# Export experimental functions
export ProcessAlgorithms, OutputStatsAlgorithms

end # module
