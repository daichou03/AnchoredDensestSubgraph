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

const DUMMY_SEED = -1

# V does not need to be relevant - in which case assign it to be something <= 0 to indicate this.
function GetGenericSeedReport(B::SparseMatrixCSC, V::Int64, R::Vector{Int64})
    inducedDS = GlobalMaximumDensity(B[R,R])
    localDS = StronglyLocalMaximumDensity(B, R, inducedDS)
    degree = V <= 0 ? -1 : GetDegree(B, V)
    rSeed(V, R, degree, GetVolume(B, R), GetInducedVolume(B, R), inducedDS.alpha_star, length(inducedDS.source_nodes), localDS.alpha_star, length(localDS.source_nodes))
end

function GetGenericSeedReportV2(B::SparseMatrixCSC, V::Int64, R::Vector{Int64})
    inducedDS = GlobalMaximumDensity(B[R,R])
    localDS = StronglyLocalMaximumDensity(B, R, inducedDS)
    degree = V <= 0 ? -1 : GetDegree(B, V)
    rSeedV2(V, R, degree, GetVolume(B, R), GetInducedVolume(B, R), densestSubgraph(inducedDS.alpha_star, R[inducedDS.source_nodes]) , localDS)
end

#----------------------------
# Densest subgraph size ratio
#----------------------------

function GetDensestSubgraphRatioSize(R::SparseMatrixCSC)
    N = size(R, 1)
    return length(GlobalMaximumDensity(R).source_nodes) / N
end
