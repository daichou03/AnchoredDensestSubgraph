using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using StatsBase # TODO: To install
using Random
using Base
include("maxflow.jl")
include("Helper_io.jl")
include("Graph_utils_yd.jl")
include("Core_algorithm_yd.jl")
include("Test_utils_yd.jl")
include("Utils.jl")

mutable struct rMinimalSeed
    R::Vector{Int64}
    inducedDS::densestSubgraph
    localDS::densestSubgraph
end

# ----------

# User input set
function GenerateUserInputSet(B::SparseMatrixCSC, V::Int64, MaxHops::Int64=2, TargetSize::Int64=8)
    pool = [V]
    for i = 1:MaxHops
        pool = GetComponentAdjacency(B, pool, true)
    end
    return GenerateUserInputSetFromPool(B, V, setdiff(pool, [V]), TargetSize)
end

function GenerateUserInputSetFromPool(B::SparseMatrixCSC, V::Int64, NeighbourPool::Vector{Int64}, TargetSize::Int64=8)
    if length(NeighbourPool) < TargetSize - 1
        pool = GetComponentAdjacency(B, NeighbourPool, true)
        return GenerateUserInputSetFromPool(B, V, setdiff(pool, [V]), TargetSize)
    end
    r = sample(NeighbourPool, TargetSize - 1, replace=false)
    append!(r, V)
    return r
end

# Reference set from user set
function GenerateReferenceSetFixedWalks(B::SparseMatrixCSC, C::Vector{Int64}, Repeats::Int64=3, Steps::Int64=2)
    r = copy(C)
    for v in C
        for i = 1:Repeats
            current = v
            for step = 1:Steps
                current = rand(GetAdjacency(B, current, false))
                r = union(r, current)
            end
        end
    end
    return r
end

