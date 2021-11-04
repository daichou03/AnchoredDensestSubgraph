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
include("CS_DBLP.jl")
include("CP_MRW.jl")
include("Test_degeneracy_yd.jl")
include("CS_generic_LA.jl")
include("CS_Evaluation_Simple.jl")

V_JW = 16028
N_JW=GetAdjacency(B,V_JW,true)
N2R_JW=GetComponentAdjacency(B,N_JW,false)
# CandidateSearch(B, V_JW, N2R_JW[1])

function CandidateSearch(B, P, V, V2)
    println(string("Degree: ", GetDegree(B, V2)))
    refinedSet = SearchNonDegRefinedSet(B, [V2], 20)
    println(string("(Not) Found refined set in # tries: ", refinedSet[2]))
    R = refinedSet[1]
    println(string("Chosen Celebrity is included in R: ", (V in R)))
    println(string("Length of R: ", length(R)))
    S_LA = LocalAnchoredDensestSubgraph(B, R).source_nodes
    S_MRW = MRW_topK(P, R, 15)
    println(string("LA: ", ReportCommunity(B, R, S_LA)))
    println(string("MRW: ", ReportCommunity(B, R, S_MRW)))
    return (R, S_LA, S_MRW)
end

function CandidateSearchForce(B, P, V, Candidates)
    for v2 in Candidates
        if (GetDegree(B, v2) < 5) || (GetDegree(B, v2) >= 15)
            continue
        end
        refinedSet = SearchNonDegRefinedSet(B, [v2], 10)
        if refinedSet[2] < 0
            continue
        end
        R = refinedSet[1]
        if !(V in R)
            continue
        end
        S_LA = LocalAnchoredDensestSubgraph(B, R).source_nodes
        S_MRW = MRW_topK(P, R, 15)
        println(string("LA: ", ReportCommunity(B, R, S_LA)))
        println(string("MRW: ", ReportCommunity(B, R, S_MRW)))
        return (v2, R, S_LA, S_MRW)
    end
    return []
end

