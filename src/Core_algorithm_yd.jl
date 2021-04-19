using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using Base
include("Memory_tracker.jl")
include("maxflow.jl") # TODO: Credit
include("Helper_io.jl")
include("Graph_utils_yd.jl")
include("Utils.jl")

# For undirected and unweighted graph.

Memory_item_GDS = "GDS"
Memory_item_GA = "GA"
Memory_item_IGA = "IGA"
Memory_item_LA = "LA"

mutable struct densestSubgraph
    alpha_star::Float64 # The minimum alpha value that can saturate all source edges
    source_nodes::Vector{Int64} # give the indices of the nodes attached to the source. Note this includes source node with index = 1, and all nodes' indices are 1 greater.
end

function GlobalDensestSubgraph(B::SparseMatrixCSC)
    stamp = RegisterFunctionStamp()

    N = size(B,1)
    # Weight for source edges
    sWeights = map(x -> GetDegree(B,x), 1:N)
    RegisterMemoryItem(Memory_item_GDS, stamp, sWeights, @varname sWeights)
    alpha_bottom = sum(sWeights) / N # Reachable
    alpha_top = maximum(sWeights) # Reachable only if alpha_bottom = alpha_top, i.e. the entire graph is a clique
    flow_alpha_minus = 0
    alpha_star = 0

    flowNetTemp = [spzeros(1,1) sparse(sWeights') spzeros(1,1);
                   spzeros(N,1) B                 sparse(repeat([alpha_bottom], N));
                   spzeros(1,N+2)]
    RegisterMemoryItem(Memory_item_GDS, stamp, flowNetTemp, @varname flowNetTemp)

    if FlowNetAlpha(flowNetTemp, alpha_bottom).flowvalue >= sum(sWeights) - 1e-6
        alpha_star = alpha_bottom
        flow_alpha_minus = FlowNetAlpha(flowNetTemp, alpha_star - 1 / (N * (N+1)))
        RegisterMemoryItem(Memory_item_GDS, stamp, flow_alpha_minus, @varname flow_alpha_minus)
    else
        while alpha_top - alpha_bottom >= 1 / (N * (N+1))
            alpha = (alpha_bottom + alpha_top) / 2
            F = FlowNetAlpha(flowNetTemp, alpha)
            RegisterMemoryItem(Memory_item_GDS, stamp, F, @varname F)
            if F.flowvalue >= sum(sWeights) - 1e-6 # YD 20201223: tolerance doesn't matter much, just don't trust alpha_top
                alpha_top = alpha
            else
                alpha_bottom = alpha
            end
        end
        DeregisterMemoryItem(Memory_item_GDS, stamp, @varname F)
        flow_alpha_minus = FlowNetAlpha(flowNetTemp, alpha_bottom)
        RegisterMemoryItem(Memory_item_GDS, stamp, flow_alpha_minus, @varname flow_alpha_minus)
        subgraph_length = length(flow_alpha_minus.source_nodes) - 1
        alpha_star = Float64(floor((alpha_bottom * subgraph_length) + 1) / subgraph_length)
    end

    ReclaimFunctionMemoryUsage(Memory_item_GDS, stamp)
    return densestSubgraph(alpha_star, PopSourceForFlowNetworkResult(flow_alpha_minus.source_nodes))
end

function FlowNetAlpha(FlowNet::SparseMatrixCSC, alpha::Float64)
    N = size(FlowNet,1) - 2
    if FlowNet[N+1, N+2] != alpha
        for i = 2:N+1
            FlowNet[i, N+2] = alpha
        end
    end
    F = maxflowYD(FlowNet)
    return F
end

# GA
# InducedDS = GlobalDensestSubgraph(B[R,R])
function GlobalAnchoredDensestSubgraph(B::SparseMatrixCSC, R::Vector{Int64}, InducedDS::densestSubgraph, ShowTrace::Bool=false)
    stamp = RegisterFunctionStamp()
    N = size(B,1)
    # Weight for source edges
    # sWeightsR = map(x -> sum(B[x,:]), R)
    density_R = InducedDS.alpha_star # Density of the densest subgraph of R
    if density_R < 1 # 20210122: This should only happen when no vertices in R connects to each other. In which case the density should be 0, and pick no vertices other than the source.
        ReclaimFunctionMemoryUsage(Memory_item_GA, stamp)
        return InducedDS
    end
    sWeightsR = map(x -> (x in R) ? GetDegree(B,x) : 0, 1:N)
    RegisterMemoryItem(Memory_item_GA, stamp, sWeightsR, @varname sWeightsR)
    alpha_bottom = density_R # Reachable (degenerate case)
    alpha_top = length(R) # Not reachable
    flow_alpha_minus = 0
    alpha_star = 0

    flowNetTemp = [spzeros(1,1) sparse(sWeightsR') spzeros(1,1);
                   spzeros(N,1) B                  sparse(repeat([alpha_bottom], N));
                   spzeros(1,N+2)]
    RegisterMemoryItem(Memory_item_GA, stamp, flowNetTemp, @varname flowNetTemp)

    if FlowNetAlphaGA(flowNetTemp, alpha_bottom).flowvalue >= sum(sWeightsR) - 1e-6
        alpha_star = alpha_bottom
        flow_alpha_minus = FlowNetAlphaGA(flowNetTemp, alpha_star - 1 / (N * (N+1)))
        RegisterMemoryItem(Memory_item_GA, stamp, flow_alpha_minus, @varname flow_alpha_minus)
    else
        while alpha_top - alpha_bottom >= 1 / (N * (N+1))
            alpha = (alpha_bottom + alpha_top) / 2
            F = FlowNetAlphaGA(flowNetTemp, alpha)
            RegisterMemoryItem(Memory_item_GA, stamp, F, @varname F)
            if F.flowvalue >= sum(sWeightsR) - 1e-6
                alpha_top = alpha
            else
                alpha_bottom = alpha
            end
            if ShowTrace
                println(string("Current alpha: ", alpha))
            end
        end
        DeregisterMemoryItem(Memory_item_GA, stamp, @varname F)
        flow_alpha_minus = FlowNetAlphaGA(flowNetTemp, alpha_bottom) 
        RegisterMemoryItem(Memory_item_GA, stamp, flow_alpha_minus, @varname flow_alpha_minus)
        subgraph_length = length(flow_alpha_minus.source_nodes) - 1
        alpha_star = Float64((floor(alpha_bottom * subgraph_length) + 1) / subgraph_length)
    end
    ReclaimFunctionMemoryUsage(Memory_item_GA, stamp)
    return densestSubgraph(alpha_star, PopSourceForFlowNetworkResult(flow_alpha_minus.source_nodes))
end

function GlobalAnchoredDensestSubgraph(B::SparseMatrixCSC, R::Vector{Int64}, ShowTrace::Bool=false)
    inducedDS = GlobalDensestSubgraph(B[R,R])
    return GlobalAnchoredDensestSubgraph(B, R, inducedDS, ShowTrace)
end

function FlowNetAlphaGA(FlowNet::SparseMatrixCSC, alpha::Float64)
    N = size(FlowNet,1) - 2
    if FlowNet[N+1, N+2] != alpha
        for i = 2:N+1
            FlowNet[i, N+2] = alpha
        end
    end
    F = maxflowYD(FlowNet)
    return F
end

# IGA

# GlobalDegree and OrderByDegreeIndices are information global to B. Pre-calculate them as below:
# GlobalDegree = map(x -> GetDegree(B,x), 1:size(B,1))
# OrderByDegreeIndices = GetOrderByDegreeGraphIndices(B)

# InducedDS = GlobalDensestSubgraph(B[R,R])
function ImprovedGlobalAnchoredDensestSubgraph(B::SparseMatrixCSC, R::Vector{Int64},
    GlobalDegree::Vector{Int64}, OrderByDegreeIndices::Array{Tuple{Int64,Int64},1}, InducedDS::densestSubgraph)
    stamp = RegisterFunctionStamp()

    N = size(B,1)
    # Weight for source edges
    density_R = InducedDS.alpha_star # Density of the densest subgraph of R
    if density_R < 1 # 20210122: This should only happen when no vertices in R connects to each other. In which case the density should be 0, and pick no vertices other than the source.
        ReclaimFunctionMemoryUsage(Memory_item_IGA, stamp)
        return InducedDS
    end
    sWeightsR = map(x -> (x in R) ? GlobalDegree[x] : 0, 1:N)
    RegisterMemoryItem(Memory_item_IGA, stamp, sWeightsR, @varname sWeightsR)
    volume_R = sum(sWeightsR)
    overdensed = GetOverdensedNodes(N, OrderByDegreeIndices, volume_R)
    RegisterMemoryItem(Memory_item_IGA, stamp, overdensed, @varname overdensed)
    rToOMatrix = B[overdensed, setdiff(1:N,overdensed)]
    RegisterMemoryItem(Memory_item_IGA, stamp, rToOMatrix, @varname rToOMatrix)
    rToOWeights = map(x -> GetDegree(rToOMatrix, x), 1:(N-length(overdensed)))
    RegisterMemoryItem(Memory_item_IGA, stamp, rToOWeights, @varname rToOWeights)
    BProp = B[setdiff(1:N,overdensed), setdiff(1:N,overdensed)]
    RegisterMemoryItem(Memory_item_IGA, stamp, BProp, @varname BProp)
    sWeightsRProp = sWeightsR[setdiff(1:N,overdensed)]
    RegisterMemoryItem(Memory_item_IGA, stamp, sWeightsRProp, @varname sWeightsRProp)

    alpha_bottom = density_R # Reachable
    alpha_top = length(R) # Not reachable
    flow_alpha_minus = 0
    alpha_star = 0

    NProp = size(BProp, 1)
    flowNetTemp = [spzeros(1,1)     sparse(sWeightsRProp') spzeros(1,1);
                   spzeros(NProp,1) BProp              sparse(repeat([alpha_bottom], NProp) + rToOWeights);
                   spzeros(1,NProp+2)]
    RegisterMemoryItem(Memory_item_IGA, stamp, flowNetTemp, @varname flowNetTemp)

    # YD: Just merge the super node with sink. Also ignore any directed edges from it to regular nodes.
    if FlowNetAlphaIGA(flowNetTemp, alpha_bottom, rToOWeights).flowvalue >= sum(sWeightsR) - 1e-6
        alpha_star = alpha_bottom
        flow_alpha_minus = FlowNetAlphaIGA(flowNetTemp, alpha_star - 1 / (N * (N+1)), rToOWeights)
        RegisterMemoryItem(Memory_item_IGA, stamp, flow_alpha_minus, @varname flow_alpha_minus)
    else
        while alpha_top - alpha_bottom >= 1 / (N * (N+1))
            alpha = (alpha_bottom + alpha_top) / 2
            F = FlowNetAlphaIGA(flowNetTemp, alpha, rToOWeights)
            RegisterMemoryItem(Memory_item_IGA, stamp, F, @varname F)
            if F.flowvalue >= sum(sWeightsR) - 1e-6
                alpha_top = alpha
            else
                alpha_bottom = alpha
            end
        end
        DeregisterMemoryItem(Memory_item_IGA, stamp, @varname F)
        flow_alpha_minus = FlowNetAlphaIGA(flowNetTemp, alpha_bottom, rToOWeights)
        RegisterMemoryItem(Memory_item_IGA, stamp, flow_alpha_minus, @varname flow_alpha_minus)
        subgraph_length = length(flow_alpha_minus.source_nodes) - 1
        alpha_star = Float64((floor(alpha_bottom * subgraph_length) + 1) / subgraph_length)
    end

    ReclaimFunctionMemoryUsage(Memory_item_IGA, stamp)
    return densestSubgraph(alpha_star, PopSourceForFlowNetworkResult(flow_alpha_minus.source_nodes))
end

function ImprovedGlobalAnchoredDensestSubgraph(B::SparseMatrixCSC, R::Vector{Int64}, GlobalDegree::Vector{Int64}, OrderByDegreeIndices::Array{Tuple{Int64,Int64},1})
    inducedDS = GlobalDensestSubgraph(B[R,R])
    return ImprovedGlobalAnchoredDensestSubgraph(B, R, GlobalDegree, OrderByDegreeIndices, inducedDS)
end

function GetOverdensedNodes(N::Int64, OrderByDegreeIndices::Array{Tuple{Int64,Int64},1}, volume_R::Union{Int64,Float64})
    overdensed_ind_low = 1
    overdensed_ind_high = N + 1
    overdensed_ind_curr = overdensed_ind_low
    while overdensed_ind_low < overdensed_ind_high
        overdensed_ind_curr = (overdensed_ind_low + overdensed_ind_high) ÷ 2
        if OrderByDegreeIndices[overdensed_ind_curr][2] >= volume_R
            overdensed_ind_high = overdensed_ind_curr
        else
            overdensed_ind_low = overdensed_ind_curr + 1
        end
    end
    overdensed = map(x->x[1], OrderByDegreeIndices[overdensed_ind_curr:N])
    # println(string("Overdensed nodes: ", length(overdensed), " / ", N))
    return overdensed
end

function FlowNetAlphaIGA(FlowNet::SparseMatrixCSC, alpha::Float64, rToOWeights::Vector{Int64})
    N = size(FlowNet,1) - 2
    if FlowNet[N+1, N+2] != alpha + rToOWeights[N]
        for i = 1:N
            FlowNet[i+1, N+2] = alpha + rToOWeights[i]
        end
    end
    F = maxflowYD(FlowNet)
    return F
end

# InducedDS = GlobalDensestSubgraph(B[R,R])
function LocalAnchoredDensestSubgraph(B::SparseMatrixCSC, R::Vector{Int64}, InducedDS::densestSubgraph, ShowTrace::Bool=false)
    stamp = RegisterFunctionStamp()

    if InducedDS.alpha_star < 1 # 20210122: This should only happen when no vertices in R connects to each other. In which case the density should be 0, and pick no vertices other than the source.
        ReclaimFunctionMemoryUsage(Memory_item_LA, stamp)
        return InducedDS
    end
    Expanded = Int64[]
    RSorted = sort(R)
    RegisterMemoryItem(Memory_item_LA, stamp, RSorted, @varname RSorted)
    Frontier = RSorted
    alpha = 0
    S = Int64[]
    SUnion = Int64[]
    L = Int64[]
    while !isempty(Frontier)
        Expanded = union(Expanded, Frontier)
        RegisterMemoryItem(Memory_item_LA, stamp, Expanded, @varname Expanded)
        L = sort(union(L, GetComponentAdjacency(B, Frontier, true))) # GetComponentAdjacency is expensive, doing it incrementally.
        RegisterMemoryItem(Memory_item_LA, stamp, L, @varname L)
        result_S = GlobalAnchoredDensestSubgraph(B[L,L], orderedSubsetIndices(L, RSorted), InducedDS)
        RegisterMemoryItem(Memory_item_LA, stamp, result_S, @varname result_S)
        alpha = result_S.alpha_star
        S = L[result_S.source_nodes]
        RegisterMemoryItem(Memory_item_LA, stamp, S, @varname S)
        if ShowTrace
            println(densestSubgraph(result_S.alpha_star, S))
        end
        SUnion = union(SUnion, S)
        RegisterMemoryItem(Memory_item_LA, stamp, SUnion, @varname SUnion)
        Frontier = setdiff(S, Expanded)
        RegisterMemoryItem(Memory_item_LA, stamp, Frontier, @varname Frontier)
    end

    ReclaimFunctionMemoryUsage(Memory_item_LA, stamp)
    return densestSubgraph(alpha, S)
end

function LocalAnchoredDensestSubgraph(B::SparseMatrixCSC, R::Vector{Int64}, ShowTrace::Bool=false)
    inducedDS = GlobalDensestSubgraph(B[R,R])
    return LocalAnchoredDensestSubgraph(B, R, inducedDS, ShowTrace)
end

# ------
# Others
# ------
