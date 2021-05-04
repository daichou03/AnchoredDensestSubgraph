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
include("CS_generic.jl")

CS_AMAZON_FOLDER = "../CaseStudy/Amazon/"

AMAZON_META_FILE = string(CS_AMAZON_FOLDER, "raw/amazon-meta.txt")
AMAZON_META_TOTAL = 548552
AMAZON_GRAPH_FILE = string(CS_AMAZON_FOLDER, "IN/Amazon0302.in")

println("Reading Amazon data...")
B = readIN(AMAZON_GRAPH_FILE)

function RetrieveProductInfoAsArray(InfoType::String="title")
    filename = AMAZON_META_FILE
    io_read = open(AMAZON_META_FILE)
    info = emptyStringArray(AMAZON_META_TOTAL)
    flag_to_next_id = true
    current_id = -1
    while !eof(io_read)
        line = readline(io_read)
        if startswith(line, "Id:")
            current_id = parse(Int64, last(split(line, " "))) + 1 # File is 0-indexed, change to 1-indexed.
            flag_to_next_id = false
        elseif !flag_to_next_id
            if startswith(line, string("  ", InfoType))
                info[current_id] = line[length(InfoType) + 5 : end]
                flag_to_next_id = true
            end
        end
    end
    close(io_read)
    return info
end

B = readIN(AMAZON_GRAPH_FILE)
allTitles = RetrieveProductInfoAsArray("title")

function GetRefinedSetAmazon(C::Vector{Int64})
    GetRefinedSet(B, C, allTitles)
end

# V = 95485
# C = GenerateUserInputSet(B,V,1,4)