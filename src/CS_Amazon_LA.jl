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
include("CS_Amazon.jl")

println("Reading Amazon product info...")
AMAZON_META_FILE = string(CS_AMAZON_FOLDER, "raw/amazon-meta.txt")
AMAZON_META_TOTAL = 548552

function RetrieveProductInfoAsArray(InfoTypes::Vector{String}=["title","group"])
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

function RetrieveTopSales(TopNum::Int64=1000)
    io_read = open(AMAZON_META_FILE)
    info = emptyStringArray(TopNum)
    current_id = -1
    current_title = ""
    current_group = ""
    while !eof(io_read)
        line = readline(io_read)
        if startswith(line, "Id:")
            current_id = parse(Int64, last(split(line, " ")))
        elseif startswith(line, string("  ", "title"))
            current_title = line[length("title") + 5 : end]                
        elseif startswith(line, string("  ", "group"))
            current_group = line[length("group") + 5 : end]
        elseif startswith(line, string("  ", "salesrank"))
            rank = parse(Int64, line[length("salesrank") + 5 : end])
            if 1 <= rank <= TopNum
                info[rank] = join([current_id, current_title, current_group], "|")
            end
        end
    end
    close(io_read)
    return info
end

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

# Can't ensure we get expansive (non-degenerate) result every time, so for case study retry until we get an expansive example.
function SearchNonDegRefinedSet(C::Vector{Int64}, MaxRetry::Int64=100)
    for i = 1:MaxRetry
        R = GenerateReferenceSetFixedWalks(B, C)
        gds = GlobalDensestSubgraph(B[R,R]).source_nodes
        lds = LocalAnchoredDensestSubgraph(B,R).source_nodes
        if length(setdiff(lds, R[gds])) > 0
            return (R, i)
        end
    end
    return (Int64[], -1)
end
