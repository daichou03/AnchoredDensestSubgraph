using SparseArrays
using MAT
using LinearAlgebra
using StatsBase
using Random
using Base
include("io.jl")


# For warming up algorithms
SAMPLE_GRAPH = sparse([1,1,1,2,2,3,3,4,2,3,4,3,4,4,5,5], [2,3,4,3,4,4,5,5,1,1,1,2,2,3,3,4], ones(Float64, 16), 5, 5) # lobster.in
SAMPLE_GRAPH_R = [1,2]
SAMPLE_GRAPH_V = 1
