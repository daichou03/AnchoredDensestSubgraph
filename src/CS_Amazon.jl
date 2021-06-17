using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using StatsBase
using Random
using Base
include("Helper_io.jl")
include("Utils.jl")
include("Graph_utils_yd.jl")

# This part of the code should be made safe to be included by both our algorithm and any competitor algorithms.

CS_AMAZON_FOLDER = "../CaseStudy/Amazon/"

AMAZON_GRAPH_FILE = string(CS_AMAZON_FOLDER, "IN/com-amazon.ungraph.in")

println("Reading Amazon data...")
B = readIN(AMAZON_GRAPH_FILE)
