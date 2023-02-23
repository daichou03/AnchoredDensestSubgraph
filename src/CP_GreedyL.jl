# Reference: this is an implementation of Greedy-L algorithm:
# Chen, J., Zaïane, O., & Goebel, R. (2009, July). Local community identification in social networks. In 2009 International Conference on Advances in Social Network Analysis and Mining (pp. 237-242). IEEE.
# Except it can start from an arbitrary node set.

using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using StatsBase
using Random
using Base
include("Utils_io.jl")
include("Utils_graph.jl")
include("Utils.jl")

# TODO: can have an increment version to improve performance.
function LinScore(G::SparseMatrixCSC, D::Vector{Int64})
    return GetVolume(G[D,D]) / length(D)
end

function OutVolume(G::SparseMatrixCSC, D::Vector{Int64})
    return GetVolume(G, D) - GetVolume(G[D,D])
end

function NumBoundaryNodes(G::SparseMatrixCSC, D::Vector{Int64})
    return length(BoundaryNodes(G, D))
end

function LexScore(G::SparseMatrixCSC, D::Vector{Int64})
    return OutVolume(G, D) / NumBoundaryNodes(G, D)
end

function LScore(G::SparseMatrixCSC, D::Vector{Int64})
    return LinScore(G,D) / LexScore(G,D)
end

function GetExDegree(G::SparseMatrixCSC, D::Vector{Int64}, V::Int64)
    return length(setdiff(GetAdjacency(G, V, false), D))
end

function BoundaryNodes(G::SparseMatrixCSC, D::Vector{Int64})
    return filter(i->GetExDegree(G, D, i) > 0, D)
end

# Faster LScore calculation when adding one vertex as per formula (6) in the paper.
function LScore_inc(G::SparseMatrixCSC, D::Vector{Int64}, V::Int64,
        LinD::Float64, OutVolumeD::Int64, NumBoundaryNodesD::Int64)
    exDeg_V = GetExDegree(G, D, V)
    inDeg_V = GetDegree(G, V) - exDeg_V
    new_numBoundary = NumBoundaryNodesD + (exDeg_V == 0 ? 0 : 1)
    for i in intersect(GetAdjacency(G, V, false), D)
        if GetExDegree(G, D, i) == 1
            new_numBoundary -= 1
        end
    end
    new_LinD = (LinD * length(D) + inDeg_V * 2) / (length(D) + 1)
    new_LexD = (OutVolumeD - inDeg_V + exDeg_V) / new_numBoundary
    return new_LinD / new_LexD
end

function LScoreCommunity(G::SparseMatrixCSC, R::Vector{Int64})
    D = copy(R)
    B = BoundaryNodes(G, D)
    LinD = LinScore(G, D)
    OutVolumeD = OutVolume(G, D)
    NumBoundaryNodesD = NumBoundaryNodes(G, D)
    Lscore = -1
    new_L = LScore(G, D)
    S = GetComponentAdjacency(G,D,false)

    # Discovery Phase
    while new_L > Lscore
        Lscore = new_L
        new_n = -1
        for n_i in S
            candLScore = LScore_inc(G, D, n_i, LinD, OutVolumeD, NumBoundaryNodesD)
            if candLScore > new_L
                new_n = n_i
                new_L = candLScore
            end
        end
        if new_L > Lscore
            newLinD = LinScore(G, union(D, new_n))
            if newLinD < LinD
                # Case 2 as per paper
                S = setdiff(S, new_n)
            else
                D = union(D, new_n)
                LinD = newLinD
                OutVolumeD = OutVolume(G, D)
                NumBoundaryNodesD = NumBoundaryNodes(G, D)
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
            # Originally does not allow equal, change it to allow equal so the output is maximal
            if !(LinScore(G, D) >= LinScore(G, setdiff(D, n_i)) && LexScore(G, D) <= LexScore(G, setdiff(D, n_i)))
                D = setdiff(D, n_i)
                changed = true
            end
        end
    end
    return (D, LScore(G, D))
end
