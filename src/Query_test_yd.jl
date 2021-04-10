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
include("Test_yd.jl")

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
    IADS_overdensed_nodes::Int64
end

# ----------

DEF_USER_MAX_HOPS = 2
DEF_USER_TARGET_SIZE = 8
DEF_ANCHOR_REPEATS = 3
DEF_AHCHOR_STEPS = 2

ALG_MASK_ADS = 1
ALG_MASK_IADS = 2
ALG_MASK_SLADS = 3
ALL_ALGORITHMS = [true, true, true]

# --------------
# User input set
# --------------

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

# ---------------------------
# Reference set from user set
# ---------------------------

# Set an upper bound for node to be taken by degree.
# c_max: the maximum degree in C
# N: #nodes of the entire graph
# We cap the degree of nodes to be taken in R by c_max * X,
# Where X is smaller if c_max is close to N, and larger if c_max is small, with a bound (min_scale, max_scale).
# Otherwise, X = log(N / c_max) / log(log_scale).
mutable struct rNodeDegreeCap
    min_scale::Float64
    log_scale::Float64
    max_scale::Float64
end

#DEFAULT_R_NODE_DEGREE_CAP = rNodeDegreeCap(1.0, 2.0, 10.0)
DEFAULT_R_NODE_DEGREE_CAP = rNodeDegreeCap(8.0, 2.0, 8.0)
NULL_R_NODE_DEGREE_CAP = rNodeDegreeCap(Inf, 1.0, Inf)

function GetRNodeDegreeCap(c_max::Int64, N::Int64, RNodeDegreeCap::rNodeDegreeCap)
    scale = min(RNodeDegreeCap.max_scale, (max(RNodeDegreeCap.min_scale, log(RNodeDegreeCap.log_scale, N / c_max))))
    return floor(Int64, scale * c_max)
end

# Baseline: fixed walks
function GenerateReferenceSetFixedWalks(B::SparseMatrixCSC, C::Vector{Int64}, Repeats::Int64=DEF_ANCHOR_REPEATS, Steps::Int64=DEF_AHCHOR_STEPS,
        RNodeDegreeCap::rNodeDegreeCap=DEFAULT_R_NODE_DEGREE_CAP)
    r = copy(C)
    rDegreeCap = GetRNodeDegreeCap(maximum(map(x->GetDegree(B,x), C)), size(B,1), RNodeDegreeCap)
    for v in C
        for i = 1:Repeats
            current = v
            for step = 1:Steps
                current = rand(GetAdjacency(B, current, false))
                if GetDegree(B, current) <= rDegreeCap
                    r = union(r, current)
                end
            end
        end
    end
    return r
end

function BulkGenerateReferenceSetFixedWalks(B::SparseMatrixCSC, user_inputs::Array{Any,1}, Repeats::Int64=DEF_ANCHOR_REPEATS, Steps::Int64=DEF_AHCHOR_STEPS,
        RNodeDegreeCap::rNodeDegreeCap=DEFAULT_R_NODE_DEGREE_CAP)
    anchors = Any[]
    for i = 1:length(user_inputs)
        r = GenerateReferenceSetFixedWalks(B, user_inputs[i], Repeats, Steps, RNodeDegreeCap)
        push!(anchors, r)
    end
    return anchors
end

# Fixed target anchor size

# Copy from Test_yd.GenerateSmallRandomWalksSet with changes.
function GenerateReferenceSetTargetSize(B::SparseMatrixCSC, C::Vector{Int64}, TargetSize::Int64, MaxStep::Int64,
        RNodeDegreeCap::rNodeDegreeCap=DEFAULT_R_NODE_DEGREE_CAP, MaxRetriesMultiplier::Int64=5, ReportTrapped::Bool=false)
    if length(C) > TargetSize
        return sample(C, TargetSize, replace=false, ordered=true)
    end    
    r = copy(C)
    rDegreeCap = GetRNodeDegreeCap(maximum(map(x->GetDegree(B,x), C)), size(B,1), RNodeDegreeCap)
    step = 0
    current = rand(C)
    retries = 0
    while length(r) < TargetSize
        if step < MaxStep
            step += 1
            current = rand(GetAdjacency(B, current, false))
            if (current in r) || (GetDegree(B, current) > rDegreeCap)
                retries += 1
                if retries >= MaxRetriesMultiplier * length(C)
                    if ReportTrapped
                        println(string("[Information] Failed to finish GenerateSmallRandomWalksSet within ", MaxStep, " hops with C = ", C, ", need to allow one more step."))
                    end
                    return GenerateSmallRandomWalksSet(B, C, TargetSize, MaxStep+1, MaxRetriesMultiplier, ReportTrapped) # Allow it to explore further if can't finish
                end
            else
                retries = 0
                append!(r, current)
                if length(r) >= TargetSize
                    return r
                end
            end
        else
            step = 0
            current = rand(C)
        end
    end
    return r
end

function BulkGenerateReferenceSetTargetSize(B::SparseMatrixCSC, user_inputs::Array{Any,1}, TargetSize::Int64, MaxStep::Int64,
        RNodeDegreeCap::rNodeDegreeCap=DEFAULT_R_NODE_DEGREE_CAP, MaxRetriesMultiplier::Int64=5)
    anchors = Any[]
    for i = 1:length(user_inputs)
        r = GenerateReferenceSetTargetSize(B, user_inputs[i], TargetSize, MaxStep, RNodeDegreeCap, MaxRetriesMultiplier)
        push!(anchors, r)
    end
    return anchors
end

# -------------
# Query Process
# -------------

# Atomic query
function ProcessGlobalAnchoredDensestSubgraph(B::SparseMatrixCSC, anchors::Array{Any,1}, inducedDS_set::Array{densestSubgraph,1})
    localDS_set = Any[]
    for i = 1:length(anchors)
        R = anchors[i]
        localDS = GlobalAnchoredDensestSubgraph(B,R,inducedDS_set[i])
        push!(localDS_set, localDS)
    end
    return localDS_set
end

function ProcessImprovedGlobalAnchoredDensestSubgraph(B::SparseMatrixCSC, anchors::Array{Any,1}, inducedDS_set::Array{densestSubgraph,1}, globalDegree::Vector{Int64},
            orderByDegreeIndices::Array{Tuple{Int64,Int64},1})
    localDS_set = Any[]
    for i = 1:length(anchors)
        R = anchors[i]
        localDS = ImprovedGlobalAnchoredDensestSubgraph(B,R,globalDegree,orderByDegreeIndices,inducedDS_set[i])
        push!(localDS_set, localDS)
    end
    return localDS_set
end

function ProcessLocalAnchoredDensestSubgraph(B::SparseMatrixCSC, anchors::Array{Any,1}, inducedDS_set::Array{densestSubgraph,1})
    localDS_set = Any[]
    for i = 1:length(anchors)
        R = anchors[i]
        localDS = LocalAnchoredDensestSubgraph(B,R,inducedDS_set[i])
        push!(localDS_set, localDS)
    end
    return localDS_set
end

# Query utils
function RetrieveDataPointsFromReport(report::rMinimalSeed)
    return RetrieveDataPointsFromReport(report.R, report.induceDS, report.localDS)
end

function RetrieveDataPointsFromReport(R::Vector{Int64}, inducedDS::densestSubgraph, localDS::densestSubgraph, IADS_overdensed_nodes::Int64=0)
    R_size = length(R)
    R_induced_DS_size = length(inducedDS.source_nodes)
    ADS_size = length(localDS.source_nodes)
    expansion_size = length(setdiff(localDS.source_nodes, R))
    return dataPoint(R_size, R_induced_DS_size, ADS_size, expansion_size, IADS_overdensed_nodes)
end

function DataPointToString(dp::dataPoint)
    return string(dp.R_size, ",", dp.R_induced_DS_size, ",", dp.ADS_size, ",", dp.expansion_size, ",", dp.IADS_overdensed_nodes)
end

function DoProcessAlgorithms(B::SparseMatrixCSC, anchors::Array{Any,1}, AlgorithmMask::Vector{Bool})
    inducedDS_set = map(r -> GlobalMaximumDensity(B[r,r]), anchors)
    globalDegree = map(x -> GetDegree(B,x), 1:size(B,1))
    orderByDegreeIndices = GetOrderByDegreeGraphIndices(B)

    performances = []
    for alg_index in 1:length(AlgorithmMask)
        if AlgorithmMask[alg_index]
            if alg_index == ALG_MASK_ADS
                append!(performances, [@timed ProcessGlobalAnchoredDensestSubgraph(B, anchors, inducedDS_set)])
            elseif alg_index == ALG_MASK_IADS
                append!(performances, [@timed ProcessImprovedGlobalAnchoredDensestSubgraph(B, anchors, inducedDS_set, globalDegree, orderByDegreeIndices)])
            else
                append!(performances, [@timed ProcessLocalAnchoredDensestSubgraph(B, anchors, inducedDS_set)])
            end
        else
            append!(performances, [[0,0,0]])
        end
    end
    return (performances, inducedDS_set, globalDegree, orderByDegreeIndices)
end

function GetSoleAlgorithmIndex(AlgorithmMask::Vector{Bool})
    for alg_index in 1:length(AlgorithmMask)
        if AlgorithmMask[alg_index]
            return alg_index
        end
    end
    error("No algorithm is chosen!")
end

function DoOutputPerformanceReports(filename::String, Tests::Int64, AlgorithmMask::Vector{Bool}, performances::Array{Any,1},
        anchors::Array{Any,1}, inducedDS_set::Array{densestSubgraph,1}, globalDegree::Vector{Int64}, orderByDegreeIndices::Array{Tuple{Int64,Int64},1})
    # Write data points to file
    mkpath("../DataPoints")
    io = open(string("../DataPoints/",filename), "w")
    N = length(globalDegree)
    for i = 1:Tests
        overdensed = 0
        if AlgorithmMask[ALG_MASK_IADS]
            # Also record #overdensed nodes for IADS
            sWeightsR = map(x -> (x in anchors[i]) ? globalDegree[x] : 0, 1:N)
            volume_R = sum(sWeightsR)
            overdensed = length(GetOverdensedNodes(N, orderByDegreeIndices, volume_R))
        end
        dataPoint = RetrieveDataPointsFromReport(anchors[i], inducedDS_set[i], performances[GetSoleAlgorithmIndex(AlgorithmMask)][1][i], overdensed)
        write(io, string(DataPointToString(dataPoint),"\n"))
    end
    close(io)
    # Write time/memory reports to file
    mkpath("../PerformanceReports")
    io = open(string("../PerformanceReports/",filename), "w")
    ret = []
    for index = 2:3
        str = ""
        for alg_index in 1:length(AlgorithmMask)
            if AlgorithmMask[alg_index]
                append!(ret, performances[alg_index][index])
                if str == ""
                    str = string(performances[alg_index][index])
                else
                    str = string(str, ",", performances[alg_index][index])
                end
            end
        end
        write(io, string(str,"\n"))
    end
    close(io)
    return ret
end

# Baseline query
function PerformQueryAllAlgorithms(B::SparseMatrixCSC, Tests::Int64, DatasetName::String,
        AlgorithmMask::Vector{Bool}=ALL_ALGORITHMS, FileNameSuffix::String="",
        MaxHops::Int64=DEF_USER_MAX_HOPS, UserTargetSize::Int64=DEF_USER_TARGET_SIZE,
        Repeats::Int64=DEF_ANCHOR_REPEATS, Steps::Int64=DEF_AHCHOR_STEPS,
        RNodeDegreeCap::rNodeDegreeCap=DEFAULT_R_NODE_DEGREE_CAP)
    user_inputs = BulkGenerateUserInputSet(B, Tests, MaxHops, UserTargetSize)
    anchors = BulkGenerateReferenceSetFixedWalks(B, user_inputs, Repeats, Steps, RNodeDegreeCap)
    (performances, inducedDS_set, globalDegree, orderByDegreeIndices) = DoProcessAlgorithms(B, anchors, AlgorithmMask)
    filename = string(DatasetName, "-", Tests, "-", MaxHops, "-", UserTargetSize, "-", Repeats, "-", Steps, FileNameSuffix)
    return DoOutputPerformanceReports(filename, Tests, AlgorithmMask, performances, anchors, inducedDS_set, globalDegree, orderByDegreeIndices)
end

# Query for fixed anchor size
function PerformQueryAllAlgorithmsAnchorSizeTest(B::SparseMatrixCSC, Tests::Int64, DatasetName::String,
        MaxHops::Int64, UserTargetSize::Int64, AnchorTargetSize::Int64, Steps::Int64,
        AlgorithmMask::Vector{Bool}=ALL_ALGORITHMS, FileNameSuffix::String="-AnchorSizeTest",
        RNodeDegreeCap::rNodeDegreeCap=DEFAULT_R_NODE_DEGREE_CAP, MaxRetriesMultiplier::Int64=5)
    user_inputs = BulkGenerateUserInputSet(B, Tests, MaxHops, UserTargetSize)
    anchors = BulkGenerateReferenceSetTargetSize(B, user_inputs, AnchorTargetSize, Steps, RNodeDegreeCap, MaxRetriesMultiplier)
    (performances, inducedDS_set, globalDegree, orderByDegreeIndices) = DoProcessAlgorithms(B, anchors, AlgorithmMask)
    filename = string(DatasetName, "-", Tests, "-", MaxHops, "-", UserTargetSize, "-", AnchorTargetSize, "-", Steps, FileNameSuffix)
    return DoOutputPerformanceReports(filename, Tests, AlgorithmMask, performances, anchors, inducedDS_set, globalDegree, orderByDegreeIndices)
end

function PerformQuerySLADSAnchorSizeTest(B::SparseMatrixCSC, Tests::Int64, DatasetName::String, MaxHops::Int64, UserTargetSize::Int64,
            AnchorTargetSize::Int64, Steps::Int64, RNodeDegreeCap::rNodeDegreeCap=DEFAULT_R_NODE_DEGREE_CAP, MaxRetriesMultiplier::Int64=5)
    return PerformQueryAllAlgorithmsAnchorSizeTest(B, Tests, DatasetName, MaxHops, UserTargetSize, AnchorTargetSize, Steps,
        [false, false, true], "-AnchorSizeTest-SLADS", RNodeDegreeCap, MaxRetriesMultiplier)
end

function PerformQueryIADSSmallAnchorSizeTest(B::SparseMatrixCSC, Tests::Int64, DatasetName::String, MaxHops::Int64, UserTargetSize::Int64,
        AnchorTargetSize::Int64, Steps::Int64, RNodeDegreeCap::rNodeDegreeCap=DEFAULT_R_NODE_DEGREE_CAP, MaxRetriesMultiplier::Int64=5)
    return PerformQueryAllAlgorithmsAnchorSizeTest(B, Tests, DatasetName, MaxHops, UserTargetSize, AnchorTargetSize, Steps,
        [true, true, false], "-AnchorSizeTest-SmallIADS", RNodeDegreeCap, MaxRetriesMultiplier)
end
# --------------
# Complete Query
# --------------

# all_test_dataset_names = ["eucore","lastfm","twitch","deezer","enron","epinion"]
# ["amazon, condmat, grqc"]
# anchor_size_test_dataset_names = ["eucore","lastfm","deezer","epinion"]
# anchor_size_test_dataset_names = ["livemocha","catster"]
# chosen_dataset_names = ["eucore","fbgov","epinion","livemocha"]

function BulkPerformQueryBaseline(dataset_names::Array{String,1}, Tests::Int64)
    for ds_name in dataset_names
        println(string("Performing Query for: ", ds_name))
        dataset = readIN(string(ds_name, ".in"))
        PerformQueryAllAlgorithms(dataset, Tests, ds_name)
    end
end

function BulkPerformQuerySLADSBaseline(dataset_names::Array{String,1}, Tests::Int64)
    for ds_name in dataset_names
        println(string("Performing Query for: ", ds_name))
        dataset = readIN(string(ds_name, ".in"))
        PerformQueryAllAlgorithms(dataset, Tests, ds_name, [false, false, true])
    end
end

function BulkPerformQueryAnchorSizeTest(dataset_names::Array{String,1}, Tests::Int64)
    for ds_name in dataset_names
        println(string("Performing Anchor Size Test Query for: ", ds_name))
        dataset = readIN(string(ds_name, ".in"))
        PerformQueryAllAlgorithmsAnchorSizeTest(dataset, Tests, ds_name, 2, 2, 8, 2)
        PerformQueryAllAlgorithmsAnchorSizeTest(dataset, Tests, ds_name, 2, 4, 16, 2)
        PerformQueryAllAlgorithmsAnchorSizeTest(dataset, Tests, ds_name, 2, 8, 32, 2)
        PerformQueryAllAlgorithmsAnchorSizeTest(dataset, Tests, ds_name, 2, 16, 64, 2)
        PerformQueryAllAlgorithmsAnchorSizeTest(dataset, Tests, ds_name, 2, 32, 128, 2)
    end
end

function BulkPerformQuerySLADSAnchorSizeTest(dataset_names::Array{String,1}, Tests::Int64)
    for ds_name in dataset_names
        println(string("Performing Anchor Size Test Query for: ", ds_name))
        dataset = readIN(string(ds_name, ".in"))
        PerformQuerySLADSAnchorSizeTest(dataset, Tests, ds_name, 2, 2, 8, 2)
        PerformQuerySLADSAnchorSizeTest(dataset, Tests, ds_name, 2, 4, 16, 2)
        PerformQuerySLADSAnchorSizeTest(dataset, Tests, ds_name, 2, 8, 32, 2)
        PerformQuerySLADSAnchorSizeTest(dataset, Tests, ds_name, 2, 16, 64, 2)
        PerformQuerySLADSAnchorSizeTest(dataset, Tests, ds_name, 2, 32, 128, 2)
    end
end

function BulkPerformQueryIADSSmallAnchorSizeTest(dataset_names::Array{String,1}, Tests::Int64)
    for ds_name in dataset_names
        println(string("Performing Small Anchor Size Test Query for: ", ds_name))
        dataset = readIN(string(ds_name, ".in"))
        # PerformQueryIADSSmallAnchorSizeTest(dataset, Tests, ds_name, 1, 2, 2, 2)
        # PerformQueryIADSSmallAnchorSizeTest(dataset, Tests, ds_name, 2, 2, 3, 2)
        PerformQueryIADSSmallAnchorSizeTest(dataset, Tests, ds_name, 2, 2, 4, 2)
        PerformQueryIADSSmallAnchorSizeTest(dataset, Tests, ds_name, 2, 2, 6, 2)
        PerformQueryIADSSmallAnchorSizeTest(dataset, Tests, ds_name, 2, 2, 8, 2)
        PerformQueryIADSSmallAnchorSizeTest(dataset, Tests, ds_name, 2, 2, 10, 2)
        PerformQueryIADSSmallAnchorSizeTest(dataset, Tests, ds_name, 2, 2, 12, 2)
    end
end

# This assumes the half edge graph data file of the original exists.
# It does not perform tests on the original graph, only the half edge ones.
function BulkPerformQueryHalfEdgeTest(dataset_names::Array{String,1}, Tests::Int64, TestOriginal::Bool=false, Iteration::Integer=5, GraphSizeThreshold::Integer=32,
        AlgorithmMask::Vector{Bool}=[false, false, true])
    for ds_name in dataset_names
        println(string("Performing Half Edge Test Query for: ", ds_name))
        if TestOriginal
            println(string("Original:"))
            filename = ds_name, ".in"
            dataset = readIN(string(ds_name, ".in"))
            PerformQueryAllAlgorithms(dataset, Tests, ds_name, AlgorithmMask)
        end
        for iter = 1:Iteration
            println(string("Iteration: ", iter))
            ds_name_half_edge = string(ds_name, "-H", iter)
            filename = ds_name_half_edge, ".in"
            dataset = readIN(string(ds_name_half_edge, ".in"))
            if size(dataset, 1) < GraphSizeThreshold
                println(string("Iteration ", iter, " size smaller than ", GraphSizeThreshold, ", stop testing for sampling from ", ds_name, "."))
                break
            else
                PerformQueryAllAlgorithms(dataset, Tests, ds_name_half_edge, AlgorithmMask)
            end
        end
    end
end

function BulkPerformDegreeCapTest(dataset_names::Array{String,1}, Tests::Int64)
    for ds_name in dataset_names
        println(string("Performing Query for: ", ds_name))
        dataset = readIN(string(ds_name, ".in"))
        caps = Array{rNodeDegreeCap,1}(undef, 5)
        caps[1] = rNodeDegreeCap(1.0, 2.0, 1.0)
        caps[2] = rNodeDegreeCap(2.0, 2.0, 2.0)
        caps[3] = rNodeDegreeCap(4.0, 2.0, 4.0)
        caps[4] = rNodeDegreeCap(8.0, 2.0, 8.0)
        caps[5] = rNodeDegreeCap(16.0, 2.0, 16.0)
        for cap in caps
            PerformQueryAllAlgorithms(dataset, Tests, ds_name, [false, false, true],
                string("-cap-",cap.min_scale,"-",cap.log_scale,"-",cap.max_scale),
                DEF_USER_MAX_HOPS, DEF_USER_TARGET_SIZE, DEF_ANCHOR_REPEATS, DEF_AHCHOR_STEPS, cap)
        end
    end
end

#,"orkut","livejournal","dblp","youtube","amazon","github","astroph","condmat","grqc","hepph","hepth","brightkite","hamster","douban","gowalla"
# IADS tests:
# ["eucore", "grqc", "hepph", "github", "livemocha", "dblp", "youtube"]

# ------
# Others
# ------

# Do SLADS but returns expanded only.
function SLADSExpansionSizeOnly(B::SparseMatrixCSC, R::Vector{Int64}, inducedDS::densestSubgraph)
    if inducedDS.alpha_star < 1 # 20210122: This should only happen when no vertices in R connects to each other. In which case the density should be 0, and pick no vertices other than the source.
        return length(GetComponentAdjacency(B, R, true))
    end
    Expanded = Int64[]
    RSorted = sort(R)
    Frontier = RSorted
    alpha = 0
    S = Int64[]
    SUnion = Int64[]
    L = Int64[]
    while !isempty(Frontier)
        Expanded = union(Expanded, Frontier)
        L = sort(union(L, GetComponentAdjacency(B, Frontier, true))) # GetComponentAdjacency is expensive, doing it incrementally.
        result_S = GlobalAnchoredDensestSubgraph(B[L,L], orderedSubsetIndices(L, RSorted), inducedDS)
        alpha = result_S.alpha_star
        S = L[result_S.source_nodes]
        SUnion = union(SUnion, S)
        Frontier = setdiff(S, Expanded)
    end
    return length(L)
end

function BulkRetrieveSLADSExpansionSize(dataset_names::Array{String,1}, Tests::Int64)
    for ds_name in dataset_names
        B = readIN(string(ds_name, ".in"))

        user_inputs = BulkGenerateUserInputSet(B, Tests)
        anchors = BulkGenerateReferenceSetFixedWalks(B, user_inputs)

        inducedDS_set = map(r -> GlobalMaximumDensity(B[r,r]), anchors)
        globalDegree = map(x -> GetDegree(B,x), 1:size(B,1))
        orderByDegreeIndices = GetOrderByDegreeGraphIndices(B)

        sizes = []
        for i in 1:Tests
            append!(sizes, SLADSExpansionSizeOnly(B,anchors[i],inducedDS_set[i]))
        end
        print(string(ds_name, ","))
        println(mean(sizes))
    end
end


# -----------
# Preparation
# -----------

# println("Loading test datasets...")
# lastfm = readIN("lastfm.in")
# eucore = readIN("eucore.in")

println("Warming up each core algorithm...")
sample_graph = sparse([1,1,1,2,2,3,3,4,2,3,4,3,4,4,5,5], [2,3,4,3,4,4,5,5,1,1,1,2,2,3,3,4], ones(Float64, 16), 5, 5) # lobster.in
GlobalMaximumDensity(sample_graph)
GlobalAnchoredDensestSubgraph(sample_graph, [1,2])
ImprovedGlobalAnchoredDensestSubgraph(sample_graph, [1,2], [3,3,4,4,2], [(5,2),(1,3),(2,3),(3,4),(4,4)])
LocalAnchoredDensestSubgraph(sample_graph, [1,2])

println("Done.")