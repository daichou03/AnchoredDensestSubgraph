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
include("CS_Simple.jl")
include("CP_MRW.jl")
include("CP_FlowSeed.jl")
include("CS_Evaluation_Simple.jl")
include("Test_degeneracy_yd.jl")

# This part of the code will run our algorithm or any reimplemented baseline algorithms.

# Simple effectiveness tests

function SampleR(B::SparseMatrixCSC, Samples::Int64=100)
    function safeSampleV(B::SparseMatrixCSC, Samples::Int64=100)
        N = size(B, 1)
        vs = fill(0, Samples)
        for i = 1:Samples
            v = rand(1:N)
            while !ConnectedComponentSizeAtLeast(B, [v], CC_SIZE_THRESHOLD)
                v = rand(1:N)
            end
            vs[i] = v
        end
        return vs
    end

    vs = safeSampleV(B, Samples)
    rs = map(v->GetStepRandomWalkFixedWalks(B, [v], 15, 4, [1.0, 1.0, 1.0, 1.0], DEFAULT_R_NODE_DEGREE_CAP), vs) # 20211107: Note has 8 times degree cap from seed node?
    rs = map(v->GetStepRandomWalkFixedWalks(B, [v], 15, 4, [1.0, 1.0, 1.0, 1.0]), vs)
    return vs, rs
end

warmed_up_LA = false

function warmupLA()
    global warmed_up_LA
    if !warmed_up_LA
        println("Warming up LA...")
        LocalAnchoredDensestSubgraph(SAMPLE_GRAPH, SAMPLE_GRAPH_R)
        warmed_up_LA = true
    end
end

function SimpleLATest(B::SparseMatrixCSC, rs)
    warmupLA()
    ss = []
    times = []
    spaces = []
    TimerReset()
    for j = 1:length(rs)
        push!(ss, LocalAnchoredDensestSubgraph(B, rs[j]).source_nodes)
        push!(times, TimerLapValue())
        push!(spaces, PopMaxMemoryUsage())
    end
    return ss, times, spaces
end

warmed_up_GL = false

function warmupGL()
    global warmed_up_GL
    if !warmed_up_GL
        println("Warming up GL...")
        LScoreCommunity(SAMPLE_GRAPH, SAMPLE_GRAPH_R)
        warmed_up_GL = true
    end
end

function SimpleGLTest(B::SparseMatrixCSC, rs)
    warmupGL()
    ss = []
    times = []
    spaces = []
    PopMaxMemoryUsage()
    TimerReset()
    for j = 1:length(rs)
        push!(ss, LScoreCommunity(B, rs[j])[1])
        push!(times, TimerLapValue())
        push!(spaces, PopMaxMemoryUsage())
    end
    return ss, times, spaces
end

warmed_up_MRW = false

function warmupMRW()
    global warmed_up_MRW
    if !warmed_up_MRW
        println("Warming up MRW...")
        MRW_topK(SAMPLE_GRAPH, SAMPLE_GRAPH_V, 2)
        warmed_up_MRW = true
    end
end

function SimpleMRWTest(P::SparseMatrixCSC, vs, rs)
    warmupMRW()
    ss = []
    times = []
    spaces = []
    PopMaxMemoryUsage()
    TimerReset()
    for j = 1:length(rs)
        push!(ss, MRW_topK(P, rs[j], 100)) # 100 should be large enough, we can always truncate the result sets before reporting.
        push!(times, TimerLapValue())
        push!(spaces, PopMaxMemoryUsage())
    end
    return ss, times, spaces
end

warmed_up_FS = false

function warmupFS()
    global warmed_up_FS
    if !warmed_up_FS
        println("Warming up FS...")
        LocalCond(SAMPLE_GRAPH, SAMPLE_GRAPH_R)
        warmed_up_FS = true
    end
end

function SimpleFSTest(B, RS, PenalityR::Float64=0.0, StrongR::Vector{Int64}=Int64[], epsilon=1.0)
    warmupFS()
    res = Any[]
    times = Float64[]
    spaces = []
    PopMaxMemoryUsage()
    TimerReset()
    for j = 1:length(RS)
        push!(res, LocalCond(B, RS[j], PenalityR, StrongR, epsilon)[1])
        push!(times, TimerLapValue())
        push!(spaces, PopMaxMemoryUsage())
    end
    return res, times, spaces
end

# Code
function BulkTestExport(RegenerateR::Bool=false)
    for dataName in SIMPLE_TEST_DATA_NAMES
        println(string("Loading ", dataName, "..."))
        B = readIN(string(dataName, ".in"))
        P = toTransitionGraph(B)
        if RegenerateR
            println("Generating R:")
            vs, rs = SampleR(B, 1000)
            ExportSimpleRs(vs, rs, dataName)
        else
            println("Importing R:")
            vs, rs = ImportSimpleRs(dataName)
        end
        # LA
        println(string("Testing LA:"))
        ss_la, times_la, spaces_la = SimpleLATest(B, rs)
        ExportSimpleResults(ss_la, times_la, spaces_la, dataName, "LA")
        ReportCommunitySimple(B, rs, ss_la, times_la, spaces_la, dataName, "LA")
        # FS
        println(string("Testing FS:"))
        ss_fs, times_fs, spaces_fs = SimpleFSTest(B, rs)
        ExportSimpleResults(ss_fs, times_fs, spaces_fs, dataName, "FS")
        ReportCommunitySimple(B, rs, ss_fs, times_fs, spaces_fs, dataName, "FS")
        # MRW
        println(string("Testing MRW:"))
        ss_mrw, times_mrw, spaces_mrw = SimpleMRWTest(P, vs, rs)
        ExportSimpleResults(ss_mrw, times_mrw, spaces_mrw, dataName, "MRW")
        ReportCommunitySimple(B, rs, ss_mrw, times_mrw, spaces_mrw, dataName, "MRW")
    end
end

# dataName = "friendster"
# B = readIN(string(dataName, ".in"))
# vs, rs = ImportSimpleRs(dataName)
# ss_la, times_la, spaces_la = SimpleLATest(B, rs)
# ExportSimpleResults(ss_la, times_la, spaces_la, dataName, "LA")
# ReportCommunitySimple(B, rs, ss_la, times_la, spaces_la, dataName, "LA")
