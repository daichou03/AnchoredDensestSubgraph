using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using StatsBase
using Random
using Base
include("../ADS/flow_network.jl")
include("../Utils_io.jl")
include("../Utils_graph.jl")
include("../ADS/algorithms.jl")
include("../Utils_warmup.jl")
include("../Utils.jl")
include("../ADS/experiments/query_tests.jl")

# From C
function GetRefinedSet(B::SparseMatrixCSC, C::Vector{Int64}, Info::Array{String, 1}, Repeats::Int64=DEF_ANCHOR_REPEATS, Steps::Int64=DEF_AHCHOR_STEPS)
    R = GenerateReferenceSetFixedWalks(B,C,Repeats,Steps)
    return GetRefinedSetFromR(B, R, Info)
end

function GetRefinedSetFromR(B::SparseMatrixCSC, R::Vector{Int64}, Info::Array{String, 1})
    refined = LocalAnchoredDensestSubgraph(B,R).source_nodes
    DisplaySubset(refined, Info)
    return refined
end

function DisplaySubset(S::Vector{Int64}, Info::Array{String, 1})
    println("Info of refined set: ")
    println("------------")
    for i in S
        println(string(i, ": ", Info[i]))
    end
    println("------------")
end
