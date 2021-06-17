using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using StatsBase
using Random
using Base
include("Helper_io.jl")
include("Graph_utils_yd.jl")
include("Utils.jl")

# TODO: can have an increment version to improve performance.
function LinScore(G::SparseMatrixCSC, D::Vector{Int64})
    return GetVolume(G[D,D]) / length(D)
end

function LexScore(G::SparseMatrixCSC, D::Vector{Int64})
    numB = length(BoundaryNodes(G, D))
    return (GetVolume(G, D) - GetVolume(G[D,D])) / numB
end

function LScore(G::SparseMatrixCSC, D::Vector{Int64})
    return LinScore(G,D) / LexScore(G,D)
end

function BoundaryNodes(G::SparseMatrixCSC, D::Vector{Int64})
    return filter(i->length(setdiff(GetAdjacency(G, i), D)) > 0, D)
end

function LScoreCommunity(G::SparseMatrixCSC, R::Vector{Int64})
    D = copy(R)
    B = BoundaryNodes(G, D)
    Lscore = -1
    new_L = LScore(G, D)
    S = GetComponentAdjacency(G,D,false)

    # Discovery Phase
    while new_L > Lscore
        Lscore = new_L
        new_n = -1
        for n_i in S
            if LScore(G, union(D, n_i)) > new_L
                new_n = n_i
                new_L = LScore(G, union(D, n_i))
            end
        end
        if new_L > Lscore
            if LinScore(G, union(D, new_n)) < LinScore(G, D)
                # Case 2 as per paper
                S = setdiff(S, new_n)
            else
                D = union(D, new_n)
                B = BoundaryNodes(G, D)
                S = GetComponentAdjacency(G,D,false)
            end
        end
    end
    # Examination Phase
    changed = true
    while changed
        changed = false
        for n_i in D
            if !(LinScore(G, D) > LinScore(G, setdiff(D, n_i)) && LexScore(G, D) < LexScore(G, setdiff(D, n_i)))
                D = setdiff(D, n_i)
                changed = true
            end
        end
    end
    return (D, LScore(G, D))
end
