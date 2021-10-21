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
include("CS_Simple.jl")
include("CP_MRW.jl")
include("CS_Evaluation.jl")
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
    rs = map(v->GetStepRandomWalkFixedWalks(B, [v], 15, 4, [1.0, 1.0, 1.0, 1.0]), vs)
    return vs, rs
end

warmed_up_LA = false

function warmupLA()
    global warmed_up_LA
    if warmed_up_LA
        println("Warming up LA...")
        LocalAnchoredDensestSubgraph(SAMPLE_GRAPH, SAMPLE_GRAPH_R)
        warmed_up_LA = true
    end
end

function SimpleLATest(B::SparseMatrixCSC, rs)
    warmupLA()
    ss = []
    times = []
    TimerReset()
    for j = 1:length(rs)
        push!(ss, LocalAnchoredDensestSubgraph(B, rs[j]).source_nodes)
        push!(times, TimerLapValue())
    end
    return ss, times
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
    TimerReset()
    for j = 1:length(rs)
        push!(ss, LScoreCommunity(B, rs[j])[1])
        push!(times, TimerLapValue())
    end
    return ss, times
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
    TimerReset()
    for j = 1:length(rs)
        push!(ss, MRW_topK(P, vs[j], length(rs[j])))
        push!(times, TimerLapValue())
    end
    return ss, times
end


# Code
DATA_GRAPH_NAMES = ["amazon", "dblp", "youtube", "skitter", "livejournal", "orkut"]

function BulkTestExport()
    for dataName in DATA_GRAPH_NAMES
        B = readIN(string(dataName, ".in"))
        P = toTransitionGraph(B)
        # Generate and import vs, rs (skip if already exported)
        vs, rs = SampleR(B, 1000)
        ExportSimpleRs(vs, rs, dataName)
        # vs, rs already exported
        # vs, rs = ImportSimpleRs(dataName)
        # LA
        ss_la, times_la = SimpleLATest(B, rs)
        ExportSimpleResults(ss_la, times_la, dataName, "LA")
        # MRW
        ss_mrw, times_mrw = SimpleMRWTest(P, vs, rs)
        ExportSimpleResults(ss_mrw, times_mrw, dataName, "MRW")
    end
end

