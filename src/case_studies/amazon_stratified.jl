using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using StatsBase
using Random
using Base
include("../ADS/flow_network.jl")
include("../utils/Utils.jl")
include("../ADS/algorithms.jl")
include("generic.jl")
include("amazon.jl")
include("../competitors/home_coded/mrw.jl")
include("../competitors/home_coded/flowseed.jl")
include("../ADS/experiments/degeneracy_tests.jl")
include("generic_la.jl")

# Old stratified test (Amazon only)

# Stratified effectiveness tests

function SampleRByDegree(Indices, Samples::Int64=100)
    rs = Any[]
    for i = 1:length(Indices)
        ind_sample = StatsBase.sample(Indices[i], Samples)
        append!(rs, 0)
        rs[i] = []
        for j in 1:Samples
            append!(rs[i], 0)
            rs[i][j] = GetStepRandomWalkFixedWalks(B, [ind_sample[j]], 18, 4, [1.0, 0.7, 0.4, 0.1])
        end
    end
    return rs
end

function StratifiedLATest(RSS)
    res = Any[]
    for i = 1:length(RSS)
        append!(res, 0)
        res[i] = []
        TimerReset()
        for j = 1:length(RSS[i])
            append!(res[i], 0)
            res[i][j] = LocalAnchoredDensestSubgraph(B, RSS[i][j]).source_nodes
        end
        println(TimerLapValue())
    end
    return res
end

# S_LA = StratifiedLATest(rss)
# println("------GL below ------")
# S_GL = StratifiedGLTest(rss)

function StratifiedGLTest(RSS)
    res = Any[]
    for i = 1:length(RSS)
        append!(res, 0)
        res[i] = []
        TimerReset()
        for j = 1:length(RSS[i])
            append!(res[i], 0)
            res[i][j] = LScoreCommunity(B, RSS[i][j])[1]
        end
        println(TimerLapValue())
    end
    return res
end

function StratifiedFSTest(RSS, PenalityR::Float64=0.0, StrongR::Vector{Int64}=Int64[], epsilon=1.0)
    res = Any[]
    for i = 1:length(RSS)
        append!(res, 0)
        res[i] = []
        TimerReset()
        for j = 1:length(RSS[i])
            append!(res[i], 0)
            res[i][j] = LocalCond(B, RSS[i][j], PenalityR, StrongR, epsilon)[1]
        end
        println(TimerLapValue())
    end
    return res
end
