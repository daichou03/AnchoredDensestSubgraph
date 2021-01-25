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

function GetGenericSeedReport(B::SparseMatrixCSC, V::Int64, R::Vector{Int64})
    inducedMD = GlobalMaximumDensity(B[R,R])
    localMD = LocalMaximumDensityV2(B, R)
    rSeed(V, R, GetDegree(B, V), GetVolume(B, R), GetInducedVolume(B, R), inducedMD.alpha_star, length(inducedMD.source_nodes)-1, localMD.alpha_star, length(localMD.source_nodes)-1)
end

#----------------------------
# Densest subgraph size ratio
#----------------------------

function GetDensestSubgraphRatioSize(R::SparseMatrixCSC)
    N = size(R, 1)
    return (length(GlobalMaximumDensity(R).source_nodes) - 1) / N
end
