# This relies on competitor codes in root/CS_Competitors. That folder is ignored in git.
# There are conflicts with our Core_algorithm_yd, so do not load that together with this code.
# Credit: FlowSeed: https://github.com/nveldt/FlowSeed

include("CS_Amazon.jl")
include("CS_Simple_flowseed.jl")

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
