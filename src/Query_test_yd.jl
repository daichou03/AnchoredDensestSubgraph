# This file contains the main functions to produce the experimental results.

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

PERFORMANCE_REPORTS_DIR = "../PerformanceReports/" # Collect_results.jl uses the same constant

mutable struct rMinimalSeed
    R::Vector{Int64}
    inducedDS::densestSubgraph
    localDS::densestSubgraph
end

mutable struct dataPoint
    R_size::Int64
    R_induced_DS_size::Int64
    ADS_size::Int64 # Size of anchored densest subgraph
    expansion_size::Int64
    IGA_overdensed_nodes::Int64
end

# ----------

ANCHOR_NODES_BASE_DIR = "../AnchorNodes/"

DEF_USER_MAX_HOPS = 2
DEF_USER_TARGET_SIZE = 8
DEF_ANCHOR_REPEATS = 3
DEF_AHCHOR_STEPS = 2

ALG_MASK_GA = 1
ALG_MASK_IGA = 2
ALG_MASK_LA = 3
ALL_ALGORITHMS = [true, true, true]
LA_ONLY = [false, false, true]

CC_SIZE_THRESHOLD = 128
USER_INPUT_DEGREE_MULTIPLIER_MAX = -1.0

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

function RandomGenerateUserInputSet(B::SparseMatrixCSC, MaxHops::Int64=DEF_USER_MAX_HOPS, TargetSize::Int64=DEF_USER_TARGET_SIZE)
    N = size(B,1)
    V = rand(1:N)
    while !ConnectedComponentSizeAtLeast(B, [V], CC_SIZE_THRESHOLD)
        V = rand(1:N)
    end
    pool = [V]
    for i = 1:MaxHops
        pool = GetComponentAdjacency(B, pool, true)
    end
    if USER_INPUT_DEGREE_MULTIPLIER_MAX > 0
        degThreshold = (sum(B)/size(B,1)/2)*USER_INPUT_DEGREE_MULTIPLIER_MAX
        poolFiltered = filter(x->GetDegree(B,x)<degThreshold, pool)
        poolSize = length(pool)
        while length(poolFiltered) < DEF_USER_TARGET_SIZE
            pool = GetComponentAdjacency(B, pool, true)
            if length(pool) <= poolSize
                throw(string("This connected component don't have enough nodes that the degree < ", degThreshold))
            end
            poolSize = length(pool)
            poolFiltered = filter(x->GetDegree(B,x)<degThreshold, pool)
        end
        pool = poolFiltered
    end
    
    return GenerateUserInputSetFromPool(B, V, setdiff(pool, [V]), TargetSize)
end

function GenerateUserInputSetFromPool(B::SparseMatrixCSC, V::Int64, NeighbourPool::Vector{Int64}, TargetSize::Int64=DEF_USER_TARGET_SIZE)
    if length(NeighbourPool) < TargetSize - 1
        pool = GetComponentAdjacency(B, NeighbourPool, true)
        return GenerateUserInputSetFromPool(B, V, setdiff(pool, [V]), TargetSize)
    end
    c = StatsBase.sample(NeighbourPool, TargetSize - 1, replace=false)
    append!(c, V)
    return c
end

function BulkGenerateUserInputSet(B::SparseMatrixCSC, Tests::Int64, MaxHops::Int64=DEF_USER_MAX_HOPS, TargetSize::Int64=DEF_USER_TARGET_SIZE)
    user_inputs = Any[]
    N = size(B,1)
    for i = 1:Tests
        c = RandomGenerateUserInputSet(B, MaxHops, TargetSize)
        push!(user_inputs, c)
    end
    return user_inputs
end

# ---------------------------
# Reference set from user set
# ---------------------------

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

# Copy from Test_degeneracy_yd.GenerateSmallRandomWalksSet with changes.
function GenerateReferenceSetTargetSize(B::SparseMatrixCSC, C::Vector{Int64}, TargetSize::Int64, MaxStep::Int64,
        RNodeDegreeCap::rNodeDegreeCap=DEFAULT_R_NODE_DEGREE_CAP, MaxRetriesMultiplier::Int64=5, ReportTrapped::Bool=false)
    if length(C) > TargetSize
        return StatsBase.sample(C, TargetSize, replace=false, ordered=true)
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
                    return GenerateReferenceSetTargetSize(B, C, TargetSize, MaxStep+1, RNodeDegreeCap, MaxRetriesMultiplier, ReportTrapped) # Allow it to explore further if can't finish
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
    result_set = Any[]
    for i = 1:length(anchors)
        R = anchors[i]
        PopMaxMemoryUsage()
        result = @timed GlobalAnchoredDensestSubgraph(B,R,inducedDS_set[i])
        result = (result[1], result[2], PopMaxMemoryUsage()) # Replace memory allocation to maximum memory usage by Memory_tracker
        push!(result_set, result)
    end
    return result_set
end

function ProcessImprovedGlobalAnchoredDensestSubgraph(B::SparseMatrixCSC, anchors::Array{Any,1}, inducedDS_set::Array{densestSubgraph,1},
    globalDegree::Vector{Int64}, orderByDegreeIndices::Array{Tuple{Int64,Int64},1})
    result_set = Any[]
    for i = 1:length(anchors)
        R = anchors[i]
        PopMaxMemoryUsage()
        result = @timed ImprovedGlobalAnchoredDensestSubgraph(B,R,globalDegree,orderByDegreeIndices,inducedDS_set[i])
        result = (result[1], result[2], PopMaxMemoryUsage()) # Replace memory allocation to maximum memory usage by Memory_tracker
        push!(result_set, result)
    end
    return result_set
end

function ProcessLocalAnchoredDensestSubgraph(B::SparseMatrixCSC, anchors::Array{Any,1}, inducedDS_set::Array{densestSubgraph,1})
    result_set = Any[]
    for i = 1:length(anchors)
        R = anchors[i]
        PopMaxMemoryUsage()
        result = @timed LocalAnchoredDensestSubgraph(B,R,inducedDS_set[i])
        result = (result[1], result[2], PopMaxMemoryUsage()) # Replace memory allocation to maximum memory usage by Memory_tracker
        push!(result_set, result)
    end
    return result_set
end

# Query utils
function RetrieveDataPointsFromReport(report::rMinimalSeed)
    return RetrieveDataPointsFromReport(report.R, report.induceDS, report.localDS)
end

function RetrieveDataPointsFromReport(R::Vector{Int64}, inducedDS::densestSubgraph, localDS::densestSubgraph, IGA_overdensed_nodes::Int64=0)
    R_size = length(R)
    R_induced_DS_size = length(inducedDS.source_nodes)
    ADS_size = length(localDS.source_nodes)
    expansion_size = length(setdiff(localDS.source_nodes, R))
    return dataPoint(R_size, R_induced_DS_size, ADS_size, expansion_size, IGA_overdensed_nodes)
end

function DataPointToString(dp::dataPoint)
    return string(dp.R_size, ",", dp.R_induced_DS_size, ",", dp.ADS_size, ",", dp.expansion_size, ",", dp.IGA_overdensed_nodes)
end

function DoProcessAlgorithms(B::SparseMatrixCSC, anchors::Array{Any,1}, AlgorithmMask::Vector{Bool})
    inducedDS_set = map(r -> GlobalDensestSubgraph(B[r,r]), anchors)
    globalDegree = map(x -> GetDegree(B,x), 1:size(B,1))
    orderByDegreeIndices = GetOrderByDegreeGraphIndices(B)

    performances = []
    for alg_index in 1:length(AlgorithmMask)
        if AlgorithmMask[alg_index]
            if alg_index == ALG_MASK_GA
                append!(performances, [ProcessGlobalAnchoredDensestSubgraph(B, anchors, inducedDS_set)])
            elseif alg_index == ALG_MASK_IGA
                append!(performances, [ProcessImprovedGlobalAnchoredDensestSubgraph(B, anchors, inducedDS_set, globalDegree, orderByDegreeIndices)])
            else
                append!(performances, [ProcessLocalAnchoredDensestSubgraph(B, anchors, inducedDS_set)])
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
    return -1
end

function DoOutputPerformanceReports(filename::String, Tests::Int64, AlgorithmMask::Vector{Bool}, performances::Array{Any,1},
        anchors::Array{Any,1}, inducedDS_set::Array{densestSubgraph,1}, globalDegree::Vector{Int64}, orderByDegreeIndices::Array{Tuple{Int64,Int64},1})
    soleAlgIndex = GetSoleAlgorithmIndex(AlgorithmMask)
    if soleAlgIndex >= 0
        # Write data points to file
        mkpath("../DataPoints")
        io_dp = open(string("../DataPoints/",filename), "w")
        N = length(globalDegree)
        for i = 1:Tests
            overdensed = 0
            if AlgorithmMask[ALG_MASK_IGA]
                # Also record #overdensed nodes for IGA
                sWeightsR = map(x -> (x in anchors[i]) ? globalDegree[x] : 0, 1:N)
                volume_R = sum(sWeightsR)
                overdensed = length(GetOverdensedNodes(N, orderByDegreeIndices, volume_R))
            end
            dataPoint = RetrieveDataPointsFromReport(anchors[i], inducedDS_set[i], performances[soleAlgIndex][i][1], overdensed)
            write(io_dp, string(DataPointToString(dataPoint),"\n"))
        end
        close(io_dp)
        # Write time/memory reports to file
        mkpath("../PerformanceReports")
        io_pr = open(string(PERFORMANCE_REPORTS_DIR,filename), "w")
        ret = []
        for perf_index = 2:3
            str = ""
            for alg_index in 1:length(AlgorithmMask)
                if AlgorithmMask[alg_index]
                    sum_perf = mapreduce(x1->x1[perf_index], +, performances[alg_index])
                    append!(ret, sum_perf)
                    if str == ""
                        str = string(sum_perf)
                    else
                        str = string(str, ",", sum_perf)
                    end
                end
            end
            write(io_pr, string(str,"\n"))
        end
        close(io_pr)
        return ret
    else
        # Skip writing data points
        println("Warning: No algorithm is chosen - this is meaningless unless for benchmark purposes.")
        return []
    end   
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
# Default only run LA.
function PerformQueryAnchorSizeTest(B::SparseMatrixCSC, Tests::Int64, DatasetName::String,
        MaxHops::Int64, UserTargetSize::Int64, AnchorTargetSize::Int64, Steps::Int64,
        AlgorithmMask::Vector{Bool}=LA_ONLY, FileNameSuffix::String="-AnchorSizeTest",
        RNodeDegreeCap::rNodeDegreeCap=DEFAULT_R_NODE_DEGREE_CAP, MaxRetriesMultiplier::Int64=5)
    user_inputs = BulkGenerateUserInputSet(B, Tests, MaxHops, UserTargetSize)
    anchors = BulkGenerateReferenceSetTargetSize(B, user_inputs, AnchorTargetSize, Steps, RNodeDegreeCap, MaxRetriesMultiplier)
    (performances, inducedDS_set, globalDegree, orderByDegreeIndices) = DoProcessAlgorithms(B, anchors, AlgorithmMask)
    filename = string(DatasetName, "-", Tests, "-", MaxHops, "-", UserTargetSize, "-", AnchorTargetSize, "-", Steps, FileNameSuffix)
    return DoOutputPerformanceReports(filename, Tests, AlgorithmMask, performances, anchors, inducedDS_set, globalDegree, orderByDegreeIndices)
end

# Query for IGA performance gain vs GA with small anchor size.
function PerformQueryIGASmallAnchorSizeTest(B::SparseMatrixCSC, Tests::Int64, DatasetName::String, MaxHops::Int64, UserTargetSize::Int64,
        AnchorTargetSize::Int64, Steps::Int64, FileNameSuffix::String="-AnchorSizeTest-SmallIGA", RNodeDegreeCap::rNodeDegreeCap=DEFAULT_R_NODE_DEGREE_CAP, MaxRetriesMultiplier::Int64=5)
    return PerformQueryAnchorSizeTest(B, Tests, DatasetName, MaxHops, UserTargetSize, AnchorTargetSize, Steps,
        [true, true, false], FileNameSuffix, RNodeDegreeCap, MaxRetriesMultiplier)
end

# Query for fixed anchor size comparison, but anchor is random.
function PerformQueryLargeAnchorSizeTest(B::SparseMatrixCSC, Tests::Int64, DatasetName::String,
        AnchorTargetSize::Int64, AlgorithmMask::Vector{Bool}=ALL_ALGORITHMS, FileNameSuffix::String="-LargeAnchorSizeTest")
    anchors = []
    for i in 1:Tests
        N = size(B, 1)
        R = StatsBase.sample(1:N, AnchorTargetSize, replace=false)
        # Efficient way to ensure R has at least one edge
        nb = StatsBase.sample(GetAdjacency(B, R[1]), 1)[1]
        if !(nb in R)
            R[2] = nb
        end
        push!(anchors, R)
    end
    (performances, inducedDS_set, globalDegree, orderByDegreeIndices) = DoProcessAlgorithms(B, anchors, AlgorithmMask)
    filename = string(join([DatasetName, Tests, AnchorTargetSize], "-"), FileNameSuffix)
    return DoOutputPerformanceReports(filename, Tests, AlgorithmMask, performances, anchors, inducedDS_set, globalDegree, orderByDegreeIndices)
end

# --------------
# Complete Query
# --------------

# All_dataset_names = ["amazon","astroph","brightkite","condmat","dblp","deezer","douban","enron","epinion","fbgov","github","gowalla","grqc","hamster","hepph","hepth","lastfm","livejournal","livemocha","orkut","youtube"]
# SmallIGA_dataset_names = ["dblp","github","grqc","hepph","livemocha","youtube"]

function BulkPerformQueryBaseline(dataset_names::Array{String,1}, Tests::Int64, AlgorithmMask::Vector{Bool}=ALL_ALGORITHMS)
    for ds_name in dataset_names
        println(string("Performing Query for: ", ds_name))
        dataset = readIN(string(ds_name, ".in"))
        PerformQueryAllAlgorithms(dataset, Tests, ds_name, AlgorithmMask)
    end
end

function BatchPerformQueryAnchorSizeTest(B::SparseMatrixCSC, ds_name::String, Tests::Int64)
    PerformQueryAnchorSizeTest(B, Tests, ds_name, 2, 2, 8, 2)
    PerformQueryAnchorSizeTest(B, Tests, ds_name, 2, 4, 16, 2)
    PerformQueryAnchorSizeTest(B, Tests, ds_name, 2, 8, 32, 2)
    PerformQueryAnchorSizeTest(B, Tests, ds_name, 2, 16, 64, 2)
    PerformQueryAnchorSizeTest(B, Tests, ds_name, 2, 32, 128, 2)
end

function BulkPerformQueryAnchorSizeTest(dataset_names::Array{String,1}, Tests::Int64)
    for ds_name in dataset_names
        println(string("Performing Anchor Size Test Query for: ", ds_name))
        dataset = readIN(string(ds_name, ".in"))
        BatchPerformQueryAnchorSizeTest(dataset, ds_name, Tests)
    end
end

function BatchPerformQueryIGASmallAnchorSizeTest(B::SparseMatrixCSC, ds_name::String, Tests::Int64)
    PerformQueryIGASmallAnchorSizeTest(B, Tests, ds_name, 2, 2, 4, 2)
    PerformQueryIGASmallAnchorSizeTest(B, Tests, ds_name, 2, 2, 6, 2)
    PerformQueryIGASmallAnchorSizeTest(B, Tests, ds_name, 2, 2, 8, 2)
    PerformQueryIGASmallAnchorSizeTest(B, Tests, ds_name, 2, 2, 10, 2)
    PerformQueryIGASmallAnchorSizeTest(B, Tests, ds_name, 2, 2, 12, 2)
end

function BulkPerformQueryIGASmallAnchorSizeTest(dataset_names::Array{String,1}, Tests::Int64)
    for ds_name in dataset_names
        println(string("Performing Small Anchor Size Test Query for: ", ds_name))
        dataset = readIN(string(ds_name, ".in"))
        BatchPerformQueryIGASmallAnchorSizeTest(dataset, ds_name, Tests)
    end
end

# This assumes the half edge graph data file of the original exists.
# It does not perform tests on the original graph, only the half edge ones.
function BulkPerformQueryHalfEdgeTest(dataset_names::Array{String,1}, Tests::Int64, TestOriginal::Bool=false, Iteration::Integer=5, GraphSizeThreshold::Integer=128,
        AlgorithmMask::Vector{Bool}=LA_ONLY)
    for ds_name in dataset_names
        println(string("Performing Half Edge Test Query for: ", ds_name))
        if TestOriginal
            println(string("Original:"))
            dataset = readIN(string(ds_name, ".in"))
            PerformQueryAllAlgorithms(dataset, Tests, ds_name, AlgorithmMask)
        end
        for iter = 1:Iteration
            println(string("Iteration: ", iter))
            ds_name_half_edge = string(ds_name, "-H", iter)
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

function BatchPerformDegreeCapTest(B::SparseMatrixCSC, ds_name::String, Tests::Int64)
    caps = Array{rNodeDegreeCap,1}(undef, 5)
    caps[1] = rNodeDegreeCap(1.0, 2.0, 1.0)
    caps[2] = rNodeDegreeCap(2.0, 2.0, 2.0)
    caps[3] = rNodeDegreeCap(4.0, 2.0, 4.0)
    caps[4] = rNodeDegreeCap(8.0, 2.0, 8.0)
    caps[5] = rNodeDegreeCap(16.0, 2.0, 16.0)
    for cap in caps
        PerformQueryAllAlgorithms(B, Tests, ds_name, LA_ONLY,
            string("-cap-",cap.min_scale,"-",cap.log_scale,"-",cap.max_scale),
            DEF_USER_MAX_HOPS, DEF_USER_TARGET_SIZE, DEF_ANCHOR_REPEATS, DEF_AHCHOR_STEPS, cap)
    end
end

function BulkPerformDegreeCapTest(dataset_names::Array{String,1}, Tests::Int64)
    for ds_name in dataset_names
        println(string("Performing Query for: ", ds_name))
        dataset = readIN(string(ds_name, ".in"))
        BatchPerformDegreeCapTest(dataset, ds_name, Tests)
    end
end

function BatchPerformQueryLargeAnchorSizeTest(B::SparseMatrixCSC, Tests::Int64, ds_name::String)
    N = size(B, 1)
    target_size = 128
    while target_size < (N / 2)
        println(string("Testing size: ", target_size))
        PerformQueryLargeAnchorSizeTest(B, Tests, ds_name, target_size)
        target_size *= 2
    end
end

function BatchPerformQueryLargeAnchorSizeTest2(B::SparseMatrixCSC, TargetSize::Int64, Tests::Int64, AlgorithmMask::Vector{Bool}=ALL_ALGORITHMS, ds_name::String="orkut")
    N = size(B, 1)
    target_size = TargetSize
    for i = 1:Tests
        println(string("Testing size: ", target_size))
        PerformQueryLargeAnchorSizeTest(B, 1, ds_name, target_size, AlgorithmMask)
        target_size += 1
    end
end

# include("Query_test_yd.jl")
# B = readIN("orkut.in")
# BatchPerformQueryLargeAnchorSizeTest2(B, 131073, 5, "orkut")


# For one data graph.
function BatchPerformAllTests(B::SparseMatrixCSC, ds_name::String, Tests::Int64, LAOnly::Bool=false)
    println("Baseline test:")
    PerformQueryAllAlgorithms(B, Tests, ds_name, LAOnly ? LA_ONLY : ALL_ALGORITHMS)
    println("Anchor size test:")
    BatchPerformQueryAnchorSizeTest(B, ds_name, Tests)
    println("Cap test:")
    BatchPerformDegreeCapTest(B, ds_name, Tests)
    println(string("Nodes expanded (standard output only): ", BatchRetrieveLAExpansionSize(B, Tests)))
    println("Half edge test - Graph size for each iteration including 0 (standard output only):")
    println(string(size(B, 1), "|", div(length(B.nzval), 2)))
    BatchPerformHalfEdgeTest(ds_name, Tests)
end

# If choose false, this loads the original graph in 0.5 ^ iter-th and to try to get the LCC of rather than existing half-edged graphs.
# This is for when the graph is very large and writing half-edge graphs are very time consuming.
function BatchPerformHalfEdgeTest(ds_name::String, Tests::Int64, LoadHalfEdge::Bool=true, GraphSizeThreshold=128)
    for iter = 1:5
        B = []
        if LoadHalfEdge
            ds_name_half_edge = string(ds_name, "-H", iter)
            B = readIN(string(ds_name_half_edge, ".in"))
        else
            B = RetrieveLargestConnectedComponent(readIN(string(ds_name, ".in"), 0.5 ^ iter))
        end
        if size(B, 1) < GraphSizeThreshold
            println(string("Iteration ", iter, " size smaller than ", GraphSizeThreshold, ", stop testing for sampling from ", ds_name, "."))
            break
        else
            PerformQueryAllAlgorithms(B, Tests, string(ds_name, "-H", iter), LA_ONLY)
            println(string(size(B, 1), "|", div(length(B.nzval), 2)))
        end
    end
end

# Just easy to find
# https://www.asciiart.eu/computers/joysticks
#             ..                          ..
#          .''..''.      .--~~--.      .``..``.
#         .:''     `----'        `----'     ``:.
#         |       .    ( * )  ( * )    .       |
#        .' ....   `.  L1/L2  R1/R2  .'     _  `.
#       : .;\  /;.  :  ( * )  ( * )  :   _ (B) _ :
#      :  :) () (:   :              :   (A) _ (D) :
#       : `:/  \:'  :       B        :     (C)   :
#        :  ''''   .'   A ( * ) D    `.         :
#       .'        '   ( * ) C ( * )    `        `.
#      .'        .''.     ( * )      .``.        `.
#    .'        .'   `. (o)      (o) .'   `.        `.
#   .'       .'      `.   1(* )2   .'      `.       `.
# .'       .'         `............'         `.       `.
# `.      .' 4 BUTTON FLIGHT  YOKE W/THROTTLE `.      .'
#   `....'   Lester                       AMC   `....'


# Copy paste below to the console.
# The idea is if something is wrong, still have the loaded data graph in the memory as B ->
# ds_name = "flickr"
# include("Query_test_yd.jl")
# using Laplacians
# using Laplacians
# B = 0
# @time B = readIN(string(ds_name, ".in"))
# BatchPerformAllTests(B, ds_name, 100, length(B.nzval)>150000000)


# Generate AnchorNodes file.

# Format:
# The first line is the data graph name.
# The second line is the number of anchor sets.
# The following lines are the anchor sets.

# Example:

# eucore
# 1000
# 1,2,3,4,5,6
# 7,8,9,10,11,12
# ......

function GenerateAnchorNodesFile(ds_name::String, OutputSubDirName::String, Tests::Int64)
    dataset = readIN(string(ds_name, ".in"))

    MaxHops = DEF_USER_MAX_HOPS
    UserTargetSize = DEF_USER_TARGET_SIZE
    Repeats = DEF_ANCHOR_REPEATS
    Steps = DEF_AHCHOR_STEPS
    RNodeDegreeCap = DEFAULT_R_NODE_DEGREE_CAP

    user_inputs = BulkGenerateUserInputSet(dataset, Tests, MaxHops, UserTargetSize)
    anchors = BulkGenerateReferenceSetFixedWalks(dataset, user_inputs, Repeats, Steps, RNodeDegreeCap)

    dir = string(ANCHOR_NODES_BASE_DIR,OutputSubDirName)
    mkpath(dir)
    io = open(string(dir, ds_name, ".anchor"), "w")
    write(io, string(ds_name, "\n"))
    write(io, string(Tests), "\n")
    for anchor in anchors
        write(io, string(join(map(string, anchor),","), "\n"))
    end
    close(io)
end

function BulkGenerateAnchorNodesFile(dataset_names::Array{String,1}, OutputSubDirName::String, Tests::Int64)
    for ds_name in dataset_names
        GenerateAnchorNodesFile(ds_name, OutputSubDirName, Tests)
    end
end

# ------
# Others
# ------

# Do LA but returns the size of the area expanded only.
function LAExpansionSizeOnly(B::SparseMatrixCSC, R::Vector{Int64}, inducedDS::densestSubgraph)
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

# TODO: Specify the anchors instead to always use the same anchor node set for this test and baseline test.
# TODO: Generate separate files for each data graph instead.
function BatchRetrieveLAExpansionSize(B::SparseMatrixCSC, Tests::Int64)
    user_inputs = BulkGenerateUserInputSet(B, Tests)
    anchors = BulkGenerateReferenceSetFixedWalks(B, user_inputs)

    inducedDS_set = map(r -> GlobalDensestSubgraph(B[r,r]), anchors)
    # globalDegree = map(x -> GetDegree(B,x), 1:size(B,1))
    # orderByDegreeIndices = GetOrderByDegreeGraphIndices(B)

    sizes = []
    for i in 1:Tests
        append!(sizes, LAExpansionSizeOnly(B,anchors[i],inducedDS_set[i]))
    end
    return StatsBase.mean(sizes)
end

# Note this rewrites the entire "exp_size/exp_size_all" file!
# TODO: Should write as separate files, for now just don't use this and manually modify that file with new data instead.
# function BulkRetrieveLAExpansionSize(dataset_names::Array{String,1}, Tests::Int64)
#     io = open(string(PERFORMANCE_REPORTS_DIR, "exp_size/exp_size_all"), "w")
#     for ds_name in dataset_names
#         B = readIN(string(ds_name, ".in"))
#         mean_size = BatchRetrieveLAExpansionSize(B, Tests)
#         write(io, string(ds_name, ",", mean_size, "\n"))
#     end
#     close(io)
# end

function RetrieveLAExpansionSize(B::SparseMatrixCSC, Tests::Int64)
    user_inputs = BulkGenerateUserInputSet(B, Tests)
    anchors = BulkGenerateReferenceSetFixedWalks(B, user_inputs)

    inducedDS_set = map(r -> GlobalDensestSubgraph(B[r,r]), anchors)
    globalDegree = map(x -> GetDegree(B,x), 1:size(B,1))
    orderByDegreeIndices = GetOrderByDegreeGraphIndices(B)

    sizes = []
    for i in 1:Tests
        append!(sizes, LAExpansionSizeOnly(B,anchors[i],inducedDS_set[i]))
    end
    return sizes
end

# -----------
# Preparation
# -----------

# println("Loading test datasets...")
# lastfm = readIN("lastfm.in")
# eucore = readIN("eucore.in")

warmed_up_LA = false

function warmupLA()
    global warmed_up_LA
    if !warmed_up_LA
        println("Warming up each core algorithm...")

        GlobalDensestSubgraph(SAMPLE_GRAPH)
        GlobalAnchoredDensestSubgraph(SAMPLE_GRAPH, SAMPLE_GRAPH_R)
        ImprovedGlobalAnchoredDensestSubgraph(SAMPLE_GRAPH, SAMPLE_GRAPH_R, [3,3,4,4,2], [(5,2),(1,3),(2,3),(3,4),(4,4)])
        LocalAnchoredDensestSubgraph(SAMPLE_GRAPH, SAMPLE_GRAPH_R)
        warmed_up_LA = true
    end
end
