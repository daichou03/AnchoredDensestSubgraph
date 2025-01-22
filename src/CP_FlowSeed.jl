# Credit: FlowSeed: https://github.com/nveldt/FlowSeed

using MAT
include("CX_FlowSeed.jl")

# PenalityR: Penalty for not including R. 1.0 is their default, but 0.0 is closer to our anchored density definition.
# StrongR: flowseed can specify nodes that MUST be included.
# epsilon: 0.1 is their default, 1.0 is closer to our anchored density definition.
function LocalCond(B::SparseMatrixCSC, R::Vector{Int64}, PenalityR::Float64=0.0, StrongR::Vector{Int64}=Int64[], epsilon=1.0, MoreStats::Bool=false)
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
    result_timed = @timed FlowSeed(B,R,epsilon,pR,RinS)
    if MoreStats
        return result_timed.value, result_timed.time
    else
        return result_timed.value
    end
end