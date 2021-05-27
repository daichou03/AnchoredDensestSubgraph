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
include("CS_generic.jl")

CS_AMAZON_FOLDER = "../CaseStudy/Amazon/"

AMAZON_META_FILE = string(CS_AMAZON_FOLDER, "raw/amazon-meta.txt")
AMAZON_META_TOTAL = 548552
AMAZON_GRAPH_FILE = string(CS_AMAZON_FOLDER, "IN/com-amazon.ungraph.in")

println("Reading Amazon data...")
B = readIN(AMAZON_GRAPH_FILE)

function RetrieveProductInfoAsArray(InfoTypes::Vector{String}=["title","group"])
    filename = AMAZON_META_FILE
    io_read = open(AMAZON_META_FILE)
    info = emptyStringArray(AMAZON_META_TOTAL)
    current_id = -1
    while !eof(io_read)
        line = readline(io_read)
        if startswith(line, "Id:")
            current_id = parse(Int64, last(split(line, " ")))
        else
            for infoType in InfoTypes
                if startswith(line, string("  ", infoType))
                    if info[current_id] != ""
                        info[current_id] = string(info[current_id], " | ")
                    end
                    info[current_id] = string(info[current_id], line[length(infoType) + 5 : end])
                end
            end
        end
    end
    close(io_read)
    return info
end

B = readIN(AMAZON_GRAPH_FILE)
allTitles = RetrieveProductInfoAsArray(["title","group"])

function DisplayAdjacent(V::Int64)
    adjs = GetAdjacency(B,V)
    for adj in adjs
        println(string(adj, ": ", allTitles[adj]))
    end
    return adjs
end

function GetRefinedSetAmazon(C::Vector{Int64}, Repeats::Int64=DEF_ANCHOR_REPEATS, Steps::Int64=DEF_AHCHOR_STEPS)
    GetRefinedSet(B, C, allTitles, Repeats, Steps)
end

function GetRAmazon(C::Vector{Int64}, Repeats::Int64=DEF_ANCHOR_REPEATS, Steps::Int64=DEF_AHCHOR_STEPS)
    R = GenerateReferenceSetFixedWalks(B, C, Repeats, Steps)
    DisplaySubset(R, allTitles)
    return R
end

function GetRefinedSetFromRAmazon(R::Vector{Int64})
    GetRefinedSetFromR(B, R, allTitles)
end

# V = 9999 # "Programming and Problem Solving With C++"

# V = 43329 # "Les Miserables"
# V = 170090 # "Crime and Punishment"
# V = 19217 # "War and Peace"
# V = 72370, 72371, 79592, 100920, 102150, 170945, 229862, 261557, 266632, 291301 # "Oh My Goddess"
# V = 230805 # Forrest Gump
# V = 325872 # Off the wall

# C = GenerateUserInputSet(B,V,1,4)