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

function GetRefinedSet(B::SparseMatrixCSC, C::Vector{Int64}, Info::Array{String, 1})
    R = GenerateReferenceSetFixedWalks(B,C)
    refined = LocalAnchoredDensestSubgraph(B,R).source_nodes
    println("Info of refined set: ")
    println("------------")
    for i in refined
        println(string(Info[i]))
    end
    println("------------")
    return refined
end