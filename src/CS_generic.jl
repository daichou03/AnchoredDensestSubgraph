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
include("Query_test_yd.jl")

# From C
function GetRefinedSet(B::SparseMatrixCSC, C::Vector{Int64}, Info::Array{String, 1}, Repeats::Int64=DEF_ANCHOR_REPEATS, Steps::Int64=DEF_AHCHOR_STEPS)
    R = GenerateReferenceSetFixedWalks(B,C,Repeats,Steps)
    return GetRefinedSetFromR(B, R, Info)
end

function GetRefinedSetFromR(B::SparseMatrixCSC, R::Vector{Int64}, Info::Array{String, 1})
    refined = LocalAnchoredDensestSubgraph(B,R).source_nodes
    println("Info of refined set: ")
    println("------------")
    for i in refined
        println(string(i, ": ", Info[i]))
    end
    println("------------")
    return refined
end