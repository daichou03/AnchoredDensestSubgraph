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

CS_DBLP_CANDIDATE_FOLDER = folderString(CS_DBLP_FOLDER, "candidate")

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

# V: JW
# Candidates: JW's 1-hop or 2-hops
# Name: test series name (folder name)
# SizeMin, SizeMax: restriction on seed node's degree
# AttemptNonDegLimit: try to get non-degenerate LA result
# ForceNonDeg: skip this seed if can't get non-degenerate LA result within AttemptNonDegLimit
# DegCap: in random walk, only take nodes with degree <= max(10, deg(seedNode) ^ 2)
function CandidateSearchStore(V, Candidates, Name, SizeMin=0, SizeMax=15, AttemptNonDegLimit=10, ForceNonDeg=true, DegCap=false)
    folder = folderString(CS_DBLP_CANDIDATE_FOLDER, Name)
    mkpath(folder)
    i = 0
    for v2 in Candidates
        i += 1
        if (GetDegree(B, v2) < SizeMin) || (GetDegree(B, v2) > SizeMax)
            continue
        end
        if DegCap
            refinedSet = SearchNonDegRefinedSet(B, [v2], max(10, GetDegree(B, v2)^2), AttemptNonDegLimit)
        else
            refinedSet = SearchNonDegRefinedSet(B, [v2], AttemptNonDegLimit)
        end
        if ForceNonDeg && (refinedSet[2] < 0)
            continue
        end
        R = refinedSet[1]
        # Do not include JW in R
        if V in R
            continue
        end
        io = open(string(folder,v2,".txt"), "w")
        write(io, string(v2, "\n"))
        S_LA = LocalAnchoredDensestSubgraph(B, R).source_nodes
        S_FS = LocalCond(B, R)[1]
        S_MRW = MRW_topK(P, R, 30) # Can truncate later
        write(io, string(join(R, ","), "\n"))
        write(io, string(join(S_LA, ","), "\n"))
        write(io, string(join(S_FS, ","), "\n"))
        write(io, string(join(S_MRW, ","), "\n"))
        write(io, string(ReportCommunity(B, R, S_LA), "\n"))
        write(io, string(ReportCommunity(B, R, S_FS), "\n"))
        write(io, string(ReportCommunity(B, R, S_MRW), "\n"))
        println(string(Name, ": Node ", v2, " passed the test and report saved. Current at #", i))
        close(io)
    end
    return []
end

# Output:
# v2, R, SS (SS in 3 lines. For R and each line of SS, comma-delimited), report for each
# Folder:
# Candidate/Name/v2/

# Execution plan
# Type 1 - RW expansion
# Choose: Find from either JW's 1-hop or 2-hops
# Do not include JW in R
# Size of seed node < 15
# Choose: Either ban or not ban high degree nodes: In RW, only include nodes with degree <= max(10, deg(V) ^ 2)

# Type 2 - Manual selection: get highest collaborated
# Retrieve this from the original (multi) graph

# Stub to call that in CS_Evaluation_Simple
function ExportGraphEditorDBLP(R, Ss, Name)
    return ExportGraphEditorForDBLP(B, R, Ss, Name, folderString(CS_DBLP_FOLDER, "Single", "GraphEditor"))
end
