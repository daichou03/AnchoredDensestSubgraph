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

# For warming up algorithms
SAMPLE_GRAPH = sparse([1,1,1,2,2,3,3,4,2,3,4,3,4,4,5,5], [2,3,4,3,4,4,5,5,1,1,1,2,2,3,3,4], ones(Float64, 16), 5, 5) # lobster.in
SAMPLE_GRAPH_R = [1,2]
SAMPLE_GRAPH_V = 1
