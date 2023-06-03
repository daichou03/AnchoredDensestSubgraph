using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using StatsBase
using Random
using Base
include("maxflow.jl")
include("Utils_io.jl")
include("Utils_graph.jl")
include("Core_algorithm_yd.jl")
include("Utils_warmup.jl")
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

FOLDER_CS_DBLP_CANDIDATE = folderString(FOLDER_CS_DBLP, "candidate")

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
function CandidateSearchStore(V::Int64, Candidates, Name::String, SizeMin=0, SizeMax=20, AttemptNonDegLimit=10, ForceNonDeg=true, DegCap=true)
    folder = folderString(FOLDER_CS_DBLP_CANDIDATE, Name)
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
            refinedSet = SearchNonDegRefinedSet(B, [v2], NULL_R_NODE_DEGREE_CAP, AttemptNonDegLimit)
        end
        if ForceNonDeg && (refinedSet[2] < 0)
            continue
        end
        R = refinedSet[1]
        # Do not include JW in R
        if V in R
            continue
        end
        Ss = ComputeSS(R)
        DoWriteSearchResults(v2, R, Ss, folder)
        println(string(Name, ": Node ", v2, " passed the test and report saved. Current at #", i))
    end
    return []
end

function ComputeSS(R::Vector{Int64})
    return [LocalAnchoredDensestSubgraph(B, R).source_nodes, LocalCond(B, R)[1], MRW_topK(P, R, 30)] # Can truncate MRW later
end

function DoWriteSearchResults(V::Int64, R, Ss, folder::String)
    io = open(string(folder, V,".txt"), "w")
    write(io, string(V, "\n"))
    write(io, string(join(R, ","), "\n"))
    for j = 1:3
        write(io, string(join(Ss[j], ","), "\n"))
    end
    for j = 1:3
        write(io, string(ReportCommunity(B, R, Ss[j]), "\n"))
    end
    close(io)
end

# Generate R by collaboration

N_JW_deg = map(v->GetDegree(B, v), N_JW)
N_JW_deg = sort(collect(zip(N_JW, N_JW_deg)), by=x->x[2], rev=true)

function GetNeighbourAndWeight(V)
    vn = GetAdjacency(B, V, false)
    weights = map(v->BW[V,v], vn)
    vnw = sort(collect(zip(vn, weights)), by=x->x[2], rev=true)
    return vnw
end

# NeighbourLimit: find at most x neighbours for each hop
# CollabThreshold: only neighbours of at least x weight for each hop
function ReferenceSetByCollab(V::Int64, NeighbourLimit::Vector{Int64}=[10, 5], CollabThreshold::Vector{Float64}=[4.0, 8.0])
    return ReferenceSetByCollabRec(V, NeighbourLimit, CollabThreshold, 1, [V])
end

function ReferenceSetByCollabRec(V::Int64, NeighbourLimit::Vector{Int64}, CollabThreshold::Vector{Float64}, Layer::Int64, R::Vector{Int64})
    degThreshold = GetDegree(B, V) ^ 2
    if Layer > length(NeighbourLimit)
        return R
    end
    count = 0
    for nw in GetNeighbourAndWeight(V)
        if count >= NeighbourLimit[Layer] || nw[2] < CollabThreshold[Layer]
            break
        elseif GetDegree(B, nw[1]) <= degThreshold # Skip if too unreachable
            union!(R, nw[1])
            R = ReferenceSetByCollabRec(nw[1], NeighbourLimit, CollabThreshold, Layer+1, R)
            count += 1
        end
    end
    return R
end

# Use PrimeSeed to make this search is trackable, and are kinda random to choose old and new seeds alike.
function CollabCandidateStore(PrimeSeed::Int64=65537, StartPos::Int64=1)
    folder = folderString(FOLDER_CS_DBLP_CANDIDATE, "collab")
    mkpath(folder)
    N = size(B, 1)
    for i = StartPos:N
        v = mod(i * PrimeSeed, N) + 1
        if GetDegree(B,v) < 5 || GetDegree(B,v) >= 25
            continue
        end
        R = ReferenceSetByCollab(v)
        if length(R) < 15
            continue
        end
        Ss = ComputeSS(R)
        DoWriteSearchResults(v, R, Ss, folder)
        println(string("Node ", v, " passed the test for collab and report saved. Current at #", i))
    end
end


###################
# Export to Gephi #
###################

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

# Export chosen candidate
function ExportCandidate(CandidateName::String, CandidateV::Int64)
    folder = folderString(FOLDER_CS_DBLP_CANDIDATE, CandidateName)
    io = open(string(folder, CandidateV, ".txt"))
    v = parse(Int64, readline(io))
    RSs = []
    for i = 1:4
        push!(RSs, map(x->parse(Int64, x), split(readline(io), ",")))
    end
    ExportGraphEditorDBLP(v, RSs[1], RSs[2:4], string(CandidateName, "-", CandidateV))
    close(io)
end


# Stub to call that in CS_Evaluation_Simple
function ExportGraphEditorDBLP(V, R, Ss, Name)
    return ExportGraphEditorForDBLP(B, V, R, Ss, Name, folderString(FOLDER_CS_DBLP, "Single", "GraphEditor"))
end
