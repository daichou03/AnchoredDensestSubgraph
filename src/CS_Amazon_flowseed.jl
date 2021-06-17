# This relies on competitor codes in root/CS_Competitors. That folder is ignored in git.
# There are conflicts with our Core_algorithm_yd, so do not load that together with this code.
# Credit: FlowSeed: https://github.com/nveldt/FlowSeed

using MAT
include("CS_Amazon.jl")
include("../CS_Competitors/FlowSeed-master/algorithms/FlowSeed-1.0.jl")
include("CS_Evaluation.jl")

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

function CSTestFS(B::SparseMatrixCSC, R::Vector{Int64}, StrongR::Vector{Int64}=Int64[], epsilon=1.0, Print::Bool=true)
    S_FS = LocalCond(B,R,StrongR,epsilon)[1]
    if Print
        println(string("S_FS = ", S_FS))
        println(ReportCommunity(B,R,S_FS))
    end
    return (S_FS, ReportCommunity(B,R,S_FS))
end


# TODO: Move this to somewhere else and still refer to here. Have a document saying it is referring external codes plus credit.