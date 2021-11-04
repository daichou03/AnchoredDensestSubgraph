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
include("CP_FlowSeed.jl")
include("Test_degeneracy_yd.jl")
include("CS_generic_LA.jl")
include("CS_Evaluation_Simple.jl")

V_JW = 16028
N_JW = GetAdjacency(B, V_JW, true)
N2R_JW = GetComponentAdjacency(B, N_JW, false)
# CandidateSearch(B, V_JW, N2R_JW[1])

function CandidateSearch(V, V2)
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
    return (R, [S_LA, [], S_MRW]) # Blank for filling S_FS
end

function CandidateSearchForce(V, Candidates)
    for v2 in Candidates
        if (GetDegree(B, v2) < 3) || (GetDegree(B, v2) > 10)
            continue
        end
        refinedSet = SearchNonDegRefinedSet(B, [v2], 10)
        if refinedSet[2] < 0
            continue
        end
        R = refinedSet[1]
        # if !(V in R)
        #     continue
        # end
        S_LA = LocalAnchoredDensestSubgraph(B, R).source_nodes
        S_FS = LocalCond(B, R)[1]
        S_MRW = MRW_topK(P, R, 30) # Can truncate later
        println(string("LA: ", ReportCommunity(B, R, S_LA)))
        println(string("FS: ", ReportCommunity(B, R, S_LA)))
        println(string("MRW: ", ReportCommunity(B, R, S_MRW)))
        return (v2, R, [S_LA, S_FS, S_MRW]) # Blank for filling S_FS
    end
    return []
end

# Execution plan
# Type 1 - RW expansion
# Find from either JW's 1-hop or 2-hops
# Do not include JW in R
# Size of seed node < 15
# Ban or not ban high degree nodes: In RW, only include nodes with degree <= max(25, deg(V) ^ 2)

# Type 2 - Manual selection: get highest collaborated
# Retrieve this from the original (multi) graph

# Stub to call that in CS_Evaluation_Simple
function ExportGraphEditorDBLP(R, Ss, Name)
    return ExportGraphEditorForDBLP(B, R, Ss, Name, folderString(CS_DBLP_FOLDER, "Single", "GraphEditor"))
end
