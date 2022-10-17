using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using StatsBase # TODO: To install
using Random
using Base
include("maxflow.jl")
include("Helper_io.jl")
include("Graph_utils_yd.jl")
include("Core_algorithm_yd.jl")
include("LP_algorithm.jl")
include("Test_utils_yd.jl")
include("Utils.jl")


# Comparing LA and Local-LP-ADS#
function ProcessAlgorithms(B::SparseMatrixCSC, anchors::Array{Any,1}, SolverMask::Vector{Bool}=[true, true])
    statsAlgorithms = []
    for solver_index in 1:length(SolverMask)
        if SolverMask[solver_index]
            result_set = Array{Any}(undef, length(anchors))
            for i = 1:length(anchors)
                R = anchors[i]
                result_set[i] = DoSolveLocalADS(solver_index, B, R, true, false)
            end
            append!(statsAlgorithms, result_set)
        else
            append!(statsAlgorithms, [])
        end
    end
    return statsAlgorithms
end


warmed_up_solver = false

function WarmUpSolvers()
    global warmed_up_solver
    if not warmed_up_solver
        for i= 1:NUM_SOLVERS
            DoSolveLocalADS(i,SAMPLE_GRAPH,SAMPLE_GRAPH_R)
        end
        warmed_up_solver = true
    end
end

