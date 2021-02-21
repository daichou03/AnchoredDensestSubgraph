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
include("Test_utils_yd.jl")
include("Utils.jl")

mutable struct rMinimalSeed
    R::Vector{Int64}
    inducedDS::densestSubgraph
    localDS::densestSubgraph
end

mutable struct dataPoint
    R_size::Int64
    R_induced_DS_size::Int64
    ADS_size::Int64
    expansion_size::Int64
end

# ----------

DEF_USER_MAX_HOPS = 2
DEF_USER_TARGET_SIZE = 8
DEF_ANCHOR_REPEATS = 3
DEF_AHCHOR_STEPS = 2

# User input set
function GenerateUserInputSet(B::SparseMatrixCSC, V::Int64, MaxHops::Int64=DEF_USER_MAX_HOPS, TargetSize::Int64=DEF_USER_TARGET_SIZE)
    pool = [V]
    for i = 1:MaxHops
        pool = GetComponentAdjacency(B, pool, true)
    end
    return GenerateUserInputSetFromPool(B, V, setdiff(pool, [V]), TargetSize)
end

function GenerateUserInputSetFromPool(B::SparseMatrixCSC, V::Int64, NeighbourPool::Vector{Int64}, TargetSize::Int64=DEF_USER_TARGET_SIZE)
    if length(NeighbourPool) < TargetSize - 1
        pool = GetComponentAdjacency(B, NeighbourPool, true)
        return GenerateUserInputSetFromPool(B, V, setdiff(pool, [V]), TargetSize)
    end
    c = sample(NeighbourPool, TargetSize - 1, replace=false)
    append!(c, V)
    return c
end

function BulkGenerateUserInputSet(B::SparseMatrixCSC, Tests::Int64, MaxHops::Int64=DEF_USER_MAX_HOPS, TargetSize::Int64=DEF_USER_TARGET_SIZE)
    user_inputs = Any[]
    N = size(B,1)
    for i = 1:Tests
        V = rand(1:N)
        c = GenerateUserInputSet(B, V, MaxHops, TargetSize)
        push!(user_inputs, c)
    end
    return user_inputs
end

# Reference set from user set
function GenerateReferenceSetFixedWalks(B::SparseMatrixCSC, C::Vector{Int64}, Repeats::Int64=DEF_ANCHOR_REPEATS, Steps::Int64=DEF_AHCHOR_STEPS)
    r = copy(C)
    for v in C
        for i = 1:Repeats
            current = v
            for step = 1:Steps
                current = rand(GetAdjacency(B, current, false))
                r = union(r, current)
            end
        end
    end
    return r
end

function BulkGenerateReferenceSetFixedWalks(B::SparseMatrixCSC, user_inputs::Array{Any,1}, Repeats::Int64=DEF_ANCHOR_REPEATS, Steps::Int64=DEF_AHCHOR_STEPS)
    anchors = Any[]
    for i = 1:length(user_inputs)
        r = GenerateReferenceSetFixedWalks(B, user_inputs[i], Repeats, Steps)
        push!(anchors, r)
    end
    return anchors
end

function ProcessQueryOnReferenceSet(B::SparseMatrixCSC, anchors::Array{Any,1})
    reports = Any[]
    for i = 1:length(anchors)
        R = anchors[i]
        inducedDS = GlobalMaximumDensity(B[R,R])
        localDS = StronglyLocalMaximumDensity(B,R,inducedDS)
        push!(reports, rMinimalSeed(R, densestSubgraph(inducedDS.alpha_star, R[inducedDS.source_nodes]), localDS))
    end
    return reports
end

function ProcessLocalMaximumDensity(B::SparseMatrixCSC, anchors::Array{Any,1}, inducedDS_set::Array{densestSubgraph,1})
    localDS_set = Any[]
    for i = 1:length(anchors)
        R = anchors[i]
        localDS = LocalMaximumDensity(B,R,inducedDS_set[i])
        push!(localDS_set, localDS)
    end
    return localDS_set
end

function ProcessImprovedLocalMaximumDensity(B::SparseMatrixCSC, anchors::Array{Any,1}, inducedDS_set::Array{densestSubgraph,1}, globalDegree)
    localDS_set = Any[]
    for i = 1:length(anchors)
        R = anchors[i]
        localDS = ImprovedLocalMaximumDensity(B,R,globalDegree,inducedDS_set[i])
        push!(localDS_set, localDS)
    end
    return localDS_set
end

function ProcessStronglyLocalMaximumDensity(B::SparseMatrixCSC, anchors::Array{Any,1}, inducedDS_set::Array{densestSubgraph,1})
    localDS_set = Any[]
    for i = 1:length(anchors)
        R = anchors[i]
        localDS = StronglyLocalMaximumDensity(B,R,inducedDS_set[i])
        push!(localDS_set, localDS)
    end
    return localDS_set
end

function RetrieveDataPointsFromReport(report::rMinimalSeed)
    return RetrieveDataPointsFromReport(report.R, report.induceDS, report.localDS)
end

function RetrieveDataPointsFromReport(R::Vector{Int64}, inducedDS::densestSubgraph, localDS::densestSubgraph)
    R_size = length(R)
    R_induced_DS_size = length(inducedDS.source_nodes)
    ADS_size = length(localDS.source_nodes)
    expansion_size = length(setdiff(localDS.source_nodes, R))
    return dataPoint(R_size, R_induced_DS_size, ADS_size, expansion_size)
end

function DataPointToString(dp::dataPoint)
    return string(dp.R_size, ",", dp.R_induced_DS_size, ",", dp.ADS_size, ",", dp.expansion_size)
end

function PerformQueryAllAlgorithms(B::SparseMatrixCSC, Tests::Int64, DatasetName::String, MaxHops::Int64=DEF_USER_MAX_HOPS, TargetSize::Int64=DEF_USER_TARGET_SIZE, Repeats::Int64=DEF_ANCHOR_REPEATS, Steps::Int64=DEF_AHCHOR_STEPS)
    user_inputs = BulkGenerateUserInputSet(B, Tests, MaxHops, TargetSize)
    anchors = BulkGenerateReferenceSetFixedWalks(B, user_inputs, Repeats, Steps)
    inducedDS_set = map(r -> GlobalMaximumDensity(B[r,r]), anchors)
    globalDegree = map(x -> GetDegree(B,x), 1:size(B,1))

    timed_local = @timed ProcessLocalMaximumDensity(B, anchors, inducedDS_set)
    timed_improved_local = @timed ProcessImprovedLocalMaximumDensity(B, anchors, inducedDS_set, globalDegree)
    timed_strongly_local = @timed ProcessStronglyLocalMaximumDensity(B, anchors, inducedDS_set)
    
    inducedDS_set = map(r -> GlobalMaximumDensity(B[r,r]), anchors)
    timed_reports = @timed ProcessQueryOnReferenceSet(B, anchors)
    # Write data points to file
    filename = string(DatasetName, "-", Tests, "-", MaxHops, "-", TargetSize, "-", Repeats, "-", Steps)
    io = open(string("../DataPoints/",filename), "w")
    for i = 1:Tests
        dataPoint = RetrieveDataPointsFromReport(anchors[i], inducedDS_set[i], timed_local[1][i])
        write(io, string(DataPointToString(dataPoint),"\n"))
    end
    close(io)
    # Write time/memory reports to file
    io = open(string("../PerformanceReports/",filename), "w")
    write(io, string(timed_local[2],",",timed_improved_local[2],",",timed_strongly_local[2],"\n"))
    write(io, string(timed_local[3],",",timed_improved_local[3],",",timed_strongly_local[3],"\n"))
    close(io)
    return (timed_local[2],timed_improved_local[2],timed_strongly_local[2],timed_local[3],timed_improved_local[3],timed_strongly_local[3])
end

# -----------------
# Preload some data
# -----------------

println("Loading test datasets...")
lastfm = RetrieveLargestConnectedComponent(readIN("lastfm.in"))
eucore = RetrieveLargestConnectedComponent(readIN("eucore.in"))

println("Warming up each core algorithm...")
lobster = readIN("lobster.in", "../Example/")
# R = GetHiveRandomWalkUntilSize(eucore, 111, 32)
# globalDegree = map(x -> GetDegree(eucore,x), 1:size(eucore,1))

GlobalMaximumDensity(lobster)
LocalMaximumDensity(lobster, [1,2])
ImprovedLocalMaximumDensity(lobster, [1,2], [3,3,4,4,2])
StronglyLocalMaximumDensity(lobster, [1,2])

println("Done.")