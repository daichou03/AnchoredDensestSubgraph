# This relies on competitor codes in root/CS_Competitors. That folder is ignored in git.
# There are conflicts with our Core_algorithm_yd, so do not load that together with this code.
# Credit: FlowSeed: https://github.com/nveldt/FlowSeed

using MAT
include("../CS_Competitors/FlowSeed-master/algorithms/FlowSeed-1.0.jl")
include("Helper_io.jl")
include("Test_utils_yd.jl")
include("Utils.jl")
include("CS_generic.jl")
include("CS_Simple.jl")

# PenalityR: Penalty for not including R. 1.0 is their default, but 0.0 is closer to our anchored density definition.
# Strong nodes: flowseed can specify nodes that MUST be included.
# epsilon: 0.1 is their default, 1.0 is closer to our anchored density definition.
function LocalCond(B::SparseMatrixCSC, R::Vector{Int64}, PenalityR::Float64=0.0, StrongR::Vector{Int64}=Int64[], epsilon=1.0)
    numR = length(R)
    RinS = zeros(numR,1)
    pR = zeros(numR,1)
    for r = 1:numR
        # Strictly penalize the exclusion of nodes that are known
        # to be in the target set
        if in(R[r],StrongR)
            RinS[r] = 1
        else
            # Place a soft penalty on excluding other seed nodes
            pR[r] = PenalityR
        end
    end
   return FlowSeed(B,R,epsilon,pR,RinS)
end

function CSTestFS(B::SparseMatrixCSC, R::Vector{Int64}, StrongR::Vector{Int64}=Int64[], epsilon=1.0, Print::Bool=true)
    S_FS = LocalCond(B,R,StrongR,epsilon)[1]
    if Print
        println(string("S_FS = ", S_FS))
        println(ReportCommunity(B,R,S_FS))
    end
    return (S_FS, ReportCommunity(B,R,S_FS))
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
    TimerReset()
    for j = 1:length(RS)
        push!(res, LocalCond(B, RS[j], PenalityR, StrongR, epsilon)[1])
        push!(times, TimerLapValue())
    end
    return res, times
end

function BulkTestExportFS()
    for dataName in SIMPLE_TEST_DATA_NAMES
        println(string("Loading ", dataName, "..."))
        B = readIN(string(dataName, ".in"))
        # Generate and import vs, rs (skip if already exported)
        # vs, rs = SampleR(B, 1000)
        # ExportSimpleRs(vs, rs, dataName)
        # vs, rs already exported
        vs, rs = ImportSimpleRs(dataName)
        # FS
        println(string("Testing FS:"))
        ss_fs, times_fs = SimpleFSTest(B, rs)
        ExportSimpleResults(ss_fs, times_fs, dataName, "FS")
    end
end