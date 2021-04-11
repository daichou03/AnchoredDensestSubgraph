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

args = split(ARGS[1], ",")
data_name = string(args[1])
tests = parse(Int64, args[2])
if length(args) == 5
    doGA = string(args[3]) == "1"
    doIGA = string(args[4]) == "1"
    doLA = string(args[5]) == "1"
    BulkPerformQueryBaseline([data_name], tests, [doGA, doIGA, doLA])
elseif length(args) == 2
    BulkPerformQueryBaseline([data_name], tests)
else
    error("Must pass 2 or 5 arguments - see Memory_usage_test.jl.")
end