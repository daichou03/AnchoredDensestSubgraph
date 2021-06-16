# This relies on competitor codes in root/CS_Competitors. That folder is ignored in git.
# There are conflicts in code so interaction to our code is kept minimal.
# Credit: FlowSeed: https://github.com/nveldt/FlowSeed

using MAT
include("CS_Amazon.jl")
include("../CS_Competitors/FlowSeed-master/algorithms/FlowSeed-1.0.jl")

# Strong nodes: flowseed can specify nodes that MUST be included.
# epsilon: 0.1 is their default, 1.0 is closer to our anchored density definition.
function LocalCond(B::SparseMatrixCSC, R::Vector{Int64}, StrongR::Vector{Int64}=Int64[], epsilon=1.0)
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
            pR[r] = 1
        end
    end
   return FlowSeed(B,R,epsilon,pR,RinS)
end

# TODO: Move this to somewhere else and still refer to here. Have a document saying it is referring external codes plus credit.