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

#----------------------------
# Densest subgraph size ratio
#----------------------------

function GetDensestSubgraphRatioSize(R::SparseMatrixCSC)
    N = size(R, 1)
    return length(GlobalDensestSubgraph(R).source_nodes) / N
end
