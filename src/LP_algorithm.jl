using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using Base
using JuMP
using HiGHS
include("Helper_io.jl")
include("Graph_utils_yd.jl")
include("Utils.jl")

function SolveLPDensestSubgraph(B::SparseMatrixCSC)
    model = Model(HiGHS.Optimizer)
    edgelist = CSCToEdgeListUndirected(B)
    n = B.n
    m = length(edgelist)
    @variable(model, x[i = 1:n] >= 0)
    @variable(model, y[i = 1:m] >= 0)
    wy = Array{Int}(undef, m)
    @constraint(model, sum(x[i] for i in 1:n) <= 1)
    @objective(model, Max, sum(y))
    for i = 1:m
        u, v = edgelist[i]
        @constraint(model, y[i] <= x[u])
        @constraint(model, y[i] <= x[v])
    end
    optimize!(model)
    return densestSubgraph(objective_value(model), findall(x->value(x)>0, x))
end

# Note that "Anchored Densest Subgraph Sharp" means ADS#, which is different from ADS (which can't be LP engineered)
function SolveLPAnchoredDensestSubgraphSharp(B::SparseMatrixCSC, R::Vector{Int64})
    model = Model(HiGHS.Optimizer)
    edgelist = CSCToEdgeListUndirected(B)
    n = B.n
    m = length(edgelist)
    @variable(model, x[i = 1:n] >= 0)
    @variable(model, y[i = 1:m] >= 0)
    wy = Array{Int}(undef, m)
    @constraint(model, sum(x[i] for i in 1:n) <= 1)
    for i = 1:m
        u, v = edgelist[i]
        if (u in R) || (v in R)
            if (u in R) && (v in R)
                wy[i] = 2
            else
                wy[i] = 1
            end
            @constraint(model, y[i] <= x[u])
            @constraint(model, y[i] <= x[v])
        else
            wy[i] = -1
            @constraint(model, y[i] >= x[u] - x[v])
            @constraint(model, y[i] >= x[v] - x[u])
        end
    end
    @objective(model, Max, sum(map(*, wy, y)))
    optimize!(model)
    return densestSubgraph(objective_value(model), findall(x->value(x)>0, x))
end

function SolveLPLocalAnchoredDensestSubgraphSharp(B::SparseMatrixCSC, R::Vector{Int64}, ShowTrace::Bool=false)
    Expanded = Int64[]
    RSorted = sort(R)
    Frontier = RSorted
    alpha = 0
    S = Int64[]
    SUnion = Int64[]
    L = Int64[]
    while !isempty(Frontier)
        Expanded = union(Expanded, Frontier)
        L = sort(union(L, GetComponentAdjacency(B, Frontier, true))) # GetComponentAdjacency is expensive, doing it incrementally.
        result_S = SolveLPAnchoredDensestSubgraphSharp(B[L,L], orderedSubsetIndices(L, RSorted))
        alpha = result_S.alpha_star
        S = L[result_S.source_nodes]
        if ShowTrace
            println(densestSubgraph(result_S.alpha_star, S))
        end
        SUnion = union(SUnion, S)
        Frontier = setdiff(S, Expanded)
    end

    return densestSubgraph(alpha, S)
end