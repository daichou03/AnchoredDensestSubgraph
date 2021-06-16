using MAT
include("../CS_Competitors/FlowSeed-master/algorithms/FlowSeed-1.0.jl")

function LocalCond(B::SparseMatrixCSC, R::Vector{Int64}, epsilon=0.1)
    RN = neighborhood(B,R,1)
    numRN = length(RN)
    RNinS = zeros(numRN,1)
    pRN = zeros(numRN,1)
    for r = 1:numRN
        # Strictly penalize the exclusion of nodes that are known
        # to be in the target set
        if in(RN[r],R)
            RNinS[r] = 1
        else
            # Place a soft penalty on excluding other seed nodes
            pRN[r] = 1
        end
    end
   return FlowSeed(B,RN,epsilon,pRN,RNinS)
end

# TODO: Move this to somewhere else and still refer to here. Have a document saying it is referring external codes plus credit.