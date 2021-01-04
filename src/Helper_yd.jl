using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using Base
include("maxflow.jl") # TODO: Credit
include("Helper_io.jl")

mutable struct densestSubgraph
    alpha_star::Float64 # The minimum alpha value that can saturate all source edges
    source_nodes::Vector{Int64} # give the indices of the nodes attached to the source. Note this includes source node with index = 1, and all nodes' indices are 1 greater.
end

function GlobalMaximumDensity(B::SparseMatrixCSC)
    N = size(B,1)
    # Weight for source edges
    sWeights = map(x -> sum(B[x,:]), 1:N)
    alpha_bottom = sum(sWeights) / N # Reachable
    alpha_top = maximum(sWeights) # Reachable only if alpha_bottom = alpha_top, i.e. the entire graph is a clique
    flow_alpha_minus = 0
    alpha_star = 0

    if FlowWithAlpha(B, alpha_bottom, sWeights).flowvalue >= sum(sWeights) - 1e-6
        alpha_star = alpha_bottom
        flow_alpha_minus = FlowWithAlpha(B, alpha_star - 1 / (N * (N+1)), sWeights)
    else
        while alpha_top - alpha_bottom >= 1 / (N * (N+1))
            alpha = (alpha_bottom + alpha_top) / 2
            F = FlowWithAlpha(B, alpha, sWeights)
            if F.flowvalue >= sum(sWeights) - 1e-6 # YD 20201223: tolerance doesn't matter much, just don't trust alpha_top
                alpha_top = alpha
            else
                alpha_bottom = alpha
            end
            # println(alpha)
        end
        flow_alpha_minus = FlowWithAlpha(B, alpha_bottom, sWeights)
        subgraph_length = length(flow_alpha_minus.source_nodes) - 1
        alpha_star = Float64(floor((alpha_bottom * subgraph_length) + 1) / subgraph_length)
    end
    return densestSubgraph(alpha_star, flow_alpha_minus.source_nodes)
end

function FlowWithAlpha(B::SparseMatrixCSC, alpha::Float64, sWeights::Vector{Float64})
    N = size(B,1)
    FlowNet = [spzeros(1,1) sparse(sWeights') spzeros(1,1);
               spzeros(N,1) B                 sparse(repeat([alpha], N));
               spzeros(1,N+2)]
    F = maxflow(FlowNet, 1, N+2)
    return F
end

# R: The INDICES of R in B.
function LocalMaximumDensity(B::SparseMatrixCSC, R::Vector{Int64})
    N = size(B,1)
    # Weight for source edges
    # sWeightsR = map(x -> sum(B[x,:]), R)
    sWeightsR = map(x -> (x in R) ? sum(B[x,:]) : Float64(0), 1:N)
    density_R = GlobalMaximumDensity(B[R,R]).alpha_star
    if density_R < 1 # Case that we don't consider.
        return GlobalMaximumDensity(B[R,R])
    end
    alpha_bottom = density_R # Reachable
    alpha_top = length(R) # Not reachable
    flow_alpha_minus = 0
    alpha_star = 0

    if FlowWithAlphaLocalDensity(B, R, alpha_bottom, sWeightsR).flowvalue >= sum(sWeightsR) - 1e-6
        alpha_star = alpha_bottom
        flow_alpha_minus = FlowWithAlphaLocalDensity(B, R, alpha_star - 1 / (N * (N+1)), sWeightsR)
    else
        while alpha_top - alpha_bottom >= 1 / (N * (N+1))
            alpha = (alpha_bottom + alpha_top) / 2
            F = FlowWithAlphaLocalDensity(B, R, alpha, sWeightsR)
            if F.flowvalue >= sum(sWeightsR) - 1e-6
                alpha_top = alpha
            else
                alpha_bottom = alpha
            end
            # println(alpha)
        end
        flow_alpha_minus = FlowWithAlphaLocalDensity(B, R, alpha_bottom, sWeightsR)
        subgraph_length = length(flow_alpha_minus.source_nodes) - 1
        alpha_star = Float64((floor(alpha_bottom * subgraph_length) + 1) / subgraph_length)
    end   
    return densestSubgraph(alpha_star, flow_alpha_minus.source_nodes)
end

function FlowWithAlphaLocalDensity(B::SparseMatrixCSC, R::Vector{Int64}, alpha::Float64, sWeightsR::Vector{Float64})
    N = size(B,1)

    FlowNet = [spzeros(1,1) sparse(sWeightsR') spzeros(1,1);
               spzeros(N,1) B                  sparse(repeat([alpha], N));
               spzeros(1,N+2)]
    F = maxflow(FlowNet, 1, N+2)
    return F
end

# Compared to LocalMaximumDensity, merge all overdensed nodes that deg(v) >= 2*vol(R) to a super node based on these nodes will never be in local densest subgraph.
function ImprovedLocalMaximumDensity(B::SparseMatrixCSC, R::Vector{Int64})
    N = size(B,1)
    # Weight for source edges
    # sWeightsR = map(x -> sum(B[x,:]), R)
    sWeights = map(x -> sum(B[x,:]), 1:N)
    sWeightsR = map(x -> (x in R) ? sWeights[x] : Float64(0), 1:N)
    volume_R = sum(sWeightsR)

    overdensed = filter(x -> x > 0, map(pair -> (pair[2] >= 2*volume_R ? pair[1] : 0), zip(1:N, sWeights)))
    rToOMatrix = B[setdiff(1:N,overdensed), overdensed]
    rToOWeights = map(x -> sum(rToOMatrix[x,:]), 1:(N-length(overdensed)))
    BProp = B[setdiff(1:N,overdensed), setdiff(1:N,overdensed)]
    sWeightsRProp = sWeightsR[setdiff(1:N,overdensed)]

    density_R = GlobalMaximumDensity(B[R,R]).alpha_star
    alpha_bottom = density_R # Reachable
    alpha_top = length(R) # Not reachable
    flow_alpha_minus = 0
    alpha_star = 0

    # YD: Just merge the super node with sink. Also ignore any directed edges from it to regular nodes.
    if FlowWithAlphaImprovedLocalDensity(BProp, R, alpha_bottom, sWeightsRProp, rToOWeights).flowvalue >= sum(sWeightsR) - 1e-6
        alpha_star = alpha_bottom
        flow_alpha_minus = FlowWithAlphaImprovedLocalDensity(BProp, R, alpha_star - 1 / (N * (N+1)), sWeightsRProp, rToOWeights)
    else     
        while alpha_top - alpha_bottom >= 1 / (N * (N+1))
            alpha = (alpha_bottom + alpha_top) / 2
            F = FlowWithAlphaImprovedLocalDensity(BProp, R, alpha, sWeightsRProp, rToOWeights)
            if F.flowvalue >= sum(sWeightsR) - 1e-6
                alpha_top = alpha
            else
                alpha_bottom = alpha
            end
            # println(alpha)
        end
        flow_alpha_minus = FlowWithAlphaImprovedLocalDensity(BProp, R, alpha_bottom, sWeightsRProp, rToOWeights)
        subgraph_length = length(flow_alpha_minus.source_nodes) - 1
        alpha_star = Float64((floor(alpha_bottom * subgraph_length) + 1) / subgraph_length)
    end

    return densestSubgraph(alpha_star, flow_alpha_minus.source_nodes)
end

function FlowWithAlphaImprovedLocalDensity(BProp::SparseMatrixCSC, R::Vector{Int64}, alpha::Float64, sWeightsR::Vector{Float64}, rToOWeights::Vector{Float64})
    NProp = size(BProp,1)

#    Supernode version
#    FlowNet = [spzeros(1,1)     sparse(sWeightsR')   spzeros(1,2);
#               spzeros(NProp,1) BProp               sparse(rToOWeights') sparse(repeat([alpha], NProp));
#               spzeros(1,1)     sparse(oToRWeights) spzeros(1,1)         infValue;
#               spzeros(1,NProp+3)]
#    F = maxflow(FlowNet, 1, NProp+3)

    # Supernode = sink version
    FlowNet = [spzeros(1,1)     sparse(sWeightsR') spzeros(1,1);
               spzeros(NProp,1) BProp              sparse(repeat([alpha], NProp) + rToOWeights);
               spzeros(1,NProp+2)]
    F = maxflow(FlowNet, 1, NProp+2)
    return F
end

# function FlowWithAlphaStronglyLocalDensity(BProp::SparseMatrixCSC, R::Vector{Int64}, alpha::Float64, sWeightsR::Vector{Float64}, rToOWeights::Vector{Float64})


# ------
# R=vec([1 2])
# ------
