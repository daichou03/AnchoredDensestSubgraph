using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using StatsBase
using Random
using Base
include("maxflow.jl")
include("Helper_io.jl")
include("Graph_utils_yd.jl")
include("Core_algorithm_yd.jl")
include("Test_utils_yd.jl")
include("Utils.jl")
include("CS_generic.jl")
include("CP_MRW.jl")
include("Test_degeneracy_yd.jl")

# Can't ensure we get expansive (non-degenerate) result every time, so for case study retry until we get an expansive example.
function SearchNonDegRefinedSet(B::SparseMatrixCSC, C::Vector{Int64}, MaxRetry::Int64=100)
    for i = 1:MaxRetry
        R = GetStepRandomWalkFixedWalks(B, C, 15, 4, [1.0, 1.0, 1.0, 1.0])
        gds = GlobalDensestSubgraph(B[R,R]).source_nodes
        lds = LocalAnchoredDensestSubgraph(B,R).source_nodes
        if length(setdiff(lds, R[gds])) > 0
            return (R, i)
        end
    end
    return (Int64[], -1)
end


function CSTest(B::SparseMatrixCSC, V::Int64, Print::Bool=true)
    C = [V]
    # Try to get an expansive result for LA for some tries, otherwise just any R
    R = SearchNonDegRefinedSet(B, C, 25)[1]
    if length(R) == 0
        R = GenerateReferenceSetFixedWalks(B, C)
    end
    S_LA = LocalAnchoredDensestSubgraph(B,R).source_nodes
    S_MRW = MRW_topK(P,V,15) # Take 15 as cluster size
    if Print
        println(string("V = ", V, " # ", AMAZON_PRODUCT_INFO[V]))
        println(string("R = ", R))
        println(string("S_LA = ", S_LA))
        println(string("S_MRW = ", S_MRW))
        println(V)
        println(length(R))
        println(ReportCommunity(B,R,S_LA))
        println(ReportCommunity(B,R,S_MRW))
    end
    return (R, S_LA, S_MRW, ReportCommunity(B,R,S_LA), ReportCommunity(B,R,S_MRW))
end
