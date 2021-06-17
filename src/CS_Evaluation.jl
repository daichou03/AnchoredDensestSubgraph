using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using StatsBase
using Random
using Base
include("maxflow.jl")
include("Helper_io.jl")
include("Graph_utils_yd.jl")
include("Core_algorithm_yd.jl")
include("Test_utils_yd.jl")
include("Utils.jl")
include("CP_GreedyL.jl")

function GetDensity(B::SparseMatrixCSC, S::Vector{Int64})
    return GetVolume(B[S,S]) / length(S)
end

function GetAnchoredDensity(B::SparseMatrixCSC, R::Vector{Int64}, S::Vector{Int64})
    S_in_R = intersect(S, R)
    S_out_R = setdiff(S, S_in_R)
    return (GetVolume(B[S,S]) - GetVolume(B, S_out_R)) / length(S)
end

function GetConductance(B::SparseMatrixCSC, S::Vector{Int64})
    return 1 - GetVolume(B[S,S]) / GetVolume(B, S)
end

function GetLocalConductance(B::SparseMatrixCSC, R::Vector{Int64}, S::Vector{Int64}, epsilon::Float64=1.0)
    O_R = GetVolume(B, R) - epsilon * GetVolume(B, setdiff(S, R)) - GetVolume(B, setdiff(R, S))
    return O_R > 0 ? (GetVolume(B, S) - GetVolume(B[S,S])) / O_R : Inf
end

function GetLScore(B::SparseMatrixCSC, S::Vector{Int64})
    return LScore(B, S)
end

function GetPropRinS(R::Vector{Int64}, S::Vector{Int64})
    return 1 - length(setdiff(R, S)) / length(R)
end

function GetPropSoutR(R::Vector{Int64}, S::Vector{Int64})
    return length(setdiff(S, R)) / length(S)
end

function ReportCommunity(B::SparseMatrixCSC, R::Vector{Int64}, S::Vector{Int64})
    return join([length(S),
        GetDensity(B,S),
        GetAnchoredDensity(B,R,S),
        GetConductance(B,S),
        GetLocalConductance(B,R,S),
        GetLScore(B,S),
        GetPropRinS(R,S),
        GetPropSoutR(R,S)], "|")
end