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
include("Test_yd.jl")
include("Query_test_yd.jl")

RSS_TEST_DIR = string(ANCHOR_NODES_BASE_DIR, "Baseline/")

# args = split(ARGS[1], ",")
# data_name = string(args[1])
# tests = parse(Int64, args[2])
# if length(args) == 5
#     doGA = string(args[3]) == "1"
#     doIGA = string(args[4]) == "1"
#     doLA = string(args[5]) == "1"
#     BulkPerformQueryBaseline([data_name], tests, [doGA, doIGA, doLA])
# elseif length(args) == 2
#     BulkPerformQueryBaseline([data_name], tests)
# else
#     error("Must pass 2 or 5 arguments - see Memory_usage_test.jl.")
# end

args = split(ARGS[1], ",")
file_name = string(args[1])
anchor_index = parse(Int64, args[2])
algorithm_index = parse(Int64, args[3])
TestRSS(file_name, anchor_index, algorithm_index)

# Reads a file for time test.

# Filename:
# The first line of file contains the file name, the following rows are R. Example:
# eucore
# 1,2,3,4,5,6
# 7,8,9,10,11,12
# ......
# RIndex:
# Which row (after 1) of R to read.
# AlgorithmIndex:
# 1-3. 0 = preparation only.
function TestRSS(Filename::String, RIndex::Int64, AlgorithmIndex::Int64)
    io_test = open(string(RSS_TEST_DIR,Filename,".anchor"))
    ds_name = readline(io_test)
    anchor_text /= ""
    for i in 1:RIndex
        anchor_text = readline(io_test)
    end
    anchor = map(x->parse(Int64, x), split(anchor_text, ","))
    close(io_test)

    dataset = readIN(string(ds_name, ".in"))
    println(string("Testing RSS, data: ", ds_name, ", run: ", RIndex, ", Algorithm Index: ", AlgorithmIndex))
    algorithmMask = [false, false, false]
    if AlgorithmIndex > 0
        algorithmMask[AlgorithmIndex] = true
    end
    anchors = Any[]
    push!(anchors, anchor)
    (performances, inducedDS_set, globalDegree, orderByDegreeIndices) = DoProcessAlgorithms(dataset, anchors, algorithmMask)
    run(`ps -p $(getpid()) -o pid,comm,vsize,rss,size`)
end
