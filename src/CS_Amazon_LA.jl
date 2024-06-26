using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using StatsBase
using Random
using Base
include("maxflow.jl")
include("Utils_io.jl")
include("Utils_graph.jl")
include("Core_algorithm_yd.jl")
include("Utils_warmup.jl")
include("Utils.jl")
include("CS_generic.jl")
include("CS_Amazon.jl")
include("CS_Evaluation_Amazon.jl")
include("Test_degeneracy_yd.jl")
include("CS_generic_LA.jl")

println("Reading Amazon product info...")
AMAZON_META_FILE = string(CS_AMAZON_FOLDER, "Raw/amazon-meta.txt")
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
                        info[current_id] = string(info[current_id], "|")
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

AMAZON_PRODUCT_INFO = RetrieveProductInfoAsArray(["title","group"])

function DisplayAdjacent(V::Int64)
    adjs = GetAdjacency(B,V)
    for adj in adjs
        println(string(adj, ": ", AMAZON_PRODUCT_INFO[adj]))
    end
    return adjs
end

function GetRefinedSetAmazon(C::Vector{Int64}, Repeats::Int64=DEF_ANCHOR_REPEATS, Steps::Int64=DEF_AHCHOR_STEPS)
    GetRefinedSet(B, C, AMAZON_PRODUCT_INFO, Repeats, Steps)
end

function GetRAmazon(C::Vector{Int64}, Repeats::Int64=DEF_ANCHOR_REPEATS, Steps::Int64=DEF_AHCHOR_STEPS)
    R = GenerateReferenceSetFixedWalks(B, C, Repeats, Steps)
    DisplaySubset(R, AMAZON_PRODUCT_INFO)
    return R
end

function GetRefinedSetFromRAmazon(R::Vector{Int64})
    GetRefinedSetFromR(B, R, AMAZON_PRODUCT_INFO)
end

# V = 9999 # "Programming and Problem Solving With C++"

# V = 43329 # "Les Miserables"
# V = 170090 # "Crime and Punishment"
# V = 19217 # "War and Peace"
# V = 72370, 72371, 79592, 100920, 102150, 170945, 229862, 261557, 266632, 291301 # "Oh My Goddess"
# V = 325872 # Off the wall

# V = 230805 # "Forrest Gump"
# R = [166635, 230805, 544238, 544634, 43592, 210313, 532993, 79574, 486222, 150053, 58415, 351059, 282741, 545996, 347495, 543154, 484979, 316455, 104265, 342233, 544139]
# S_LA = [58415, 166635, 230805, 278222, 295970, 347495, 351059, 461820, 484359, 543154, 544238, 544634]
# S_FS = [166635, 230805, 544238, 544634, 43592, 210313, 532993, 79574, 486222, 150053, 58415, 351059, 282741, 545996, 347495, 543154, 484979, 316455, 104265, 342233, 544139, 49516, 94062, 278222, 295970, 461820, 468022, 484359, 485554, 512294, 543904]
# S_GD = [484979, 316455, 359454, 343531]

# V = 72370 # "Oh My Goddess"
# R = [72370, 94242, 94246, 134143, 373719, 389896, 413692, 543801, 543802, 35050, 506867, 335249, 72371, 434276, 100920, 367186, 361337, 112723, 545069, 466130, 218147, 44918, 283142, 509960, 405823, 68838, 218774, 29740, 80600, 29738, 285038]
# S_LA = [35050, 44918, 72370, 72371, 79592, 94242, 94246, 100920, 112723, 134143, 218147, 283142, 361337, 367186, 373719, 389896, 405823, 413692, 434276, 466130, 506867, 509960, 516373, 545069]
# S_FS = [72370, 94242, 94246, 134143, 373719, 389896, 413692, 543801, 543802, 35050, 506867, 335249, 72371, 434276, 100920, 367186, 361337, 112723, 545069, 466130, 218147, 44918, 283142, 509960, 405823, 68838, 218774, 29740, 80600, 29738, 285038, 79592, 134141, 219502, 245029, 262634, 289096, 344236, 428747, 516373, 544681, 547542]

# C = GenerateUserInputSet(B,V,1,4)