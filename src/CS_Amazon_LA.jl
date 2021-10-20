using Base: Bool
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
include("CP_MRW.jl")
include("CS_Evaluation.jl")
include("Test_degeneracy_yd.jl")

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
        R = GetStepRandomWalkFixedWalks(B, C, 15, 4, [1.0, 1.0, 1.0, 1.0])
        gds = GlobalDensestSubgraph(B[R,R]).source_nodes
        lds = LocalAnchoredDensestSubgraph(B,R).source_nodes
        if length(setdiff(lds, R[gds])) > 0
            return (R, i)
        end
    end
    return (Int64[], -1)
end

function CSTest(V::Int64, Print::Bool=true)
    C = GetAdjacency(B, V)
    # Try to get an expansive result for LA for some tries, otherwise just any R
    R = SearchNonDegRefinedSet(C, 25)[1]
    if length(R) == 0
        R = GenerateReferenceSetFixedWalks(B, C)
    end
    S_LA = LocalAnchoredDensestSubgraph(B,R).source_nodes
    S_MRW = MRW_topK(P,V,length(R)) # Take |R| as cluster size
    if Print
        println(string("V = ", V, " # ", allTitles[V]))
        println(string("R = ", R))
        println(string("S_LA = ", S_LA))
        println(string("S_MRW = ", S_MRW))
        println(V)
        println(length(R))
        println(ReportCommunity(B,R,S_LA))
        println(ReportCommunity(B,R,S_MRW))
    end
    return (R, S_LA, S_MRW, ReportCommunity(B,R,S_LA), ReportCommunity(B,R,S_MRW))
end

# Simple effectiveness tests

function safeSampleV(Samples::Int64=100)
    N = size(B, 1)
    vs = fill(0, Samples)
    for i = 1:Samples
        v = rand(1:N)
        while !ConnectedComponentSizeAtLeast(B, [v], CC_SIZE_THRESHOLD)
            v = rand(1:N)
        end
        vs[i] = v
    end
    return vs
end

function SampleR(Samples::Int64=100)
    vs = safeSampleV(Samples)
    rs = map(v->GetStepRandomWalkFixedWalks(B, [v], 15, 4, [1.0, 1.0, 1.0, 1.0]), vs)
    return vs, rs
end

function SimpleLATest(rs)
    res = []
    TimerReset()
    for j = 1:length(rs)
        push!(res, LocalAnchoredDensestSubgraph(B, rs[j]).source_nodes)
    end
    println(TimerLapValue())
    return res
end

warmed_up_GL = false

function warmupGL()
    global warmed_up_GL
    if !warmed_up_GL
        println("Warming up GL...")
        LScoreCommunity(SAMPLE_GRAPH, SAMPLE_GRAPH_R)
        warmed_up_GL = true

function SimpleGLTest(rs)
    warmupGL()
    res = []
    TimerReset()
    for j = 1:length(rs)
        push!(res, LScoreCommunity(B, rs[j])[1])
    end
    println(TimerLapValue())
    return res
end

warmed_up_MRW = false

function warmupMRW()
    global warmed_up_MRW
    if !warmed_up_MRW
        println("Warming up MRW...")
        MRW_topK(SAMPLE_GRAPH, SAMPLE_GRAPH_V, 2)
        warmed_up_MRW = true

function SimpleMRWTest(vs, rs)
    warmupMRW()
    res = []
    TimerReset()
    for j = 1:length(rs)
        push!(res, MRW_topK(P, vs[j], length(rs[j])))
    end
    println(TimerLapValue())
    return res
end

# TODO: Perform tests.

# Stratified effectiveness tests

function SampleRByDegree(Indices, Samples::Int64=100)
    rs = Any[]
    for i = 1:length(Indices)
        ind_sample = StatsBase.sample(Indices[i], Samples)
        append!(rs, 0)
        rs[i] = []
        for j in 1:Samples
            append!(rs[i], 0)
            rs[i][j] = GetStepRandomWalkFixedWalks(B, [ind_sample[j]], 18, 4, [1.0, 0.7, 0.4, 0.1])
        end
    end
    return rs
end

function StratifiedLATest(RSS)
    res = Any[]
    for i = 1:length(RSS)
        append!(res, 0)
        res[i] = []
        TimerReset()
        for j = 1:length(RSS[i])
            append!(res[i], 0)
            res[i][j] = LocalAnchoredDensestSubgraph(B, RSS[i][j]).source_nodes
        end
        println(TimerLapValue())
    end
    return res
end

# S_LA = StratifiedLATest(rss)
# println("------GL below ------")
# S_GL = StratifiedGLTest(rss)

function StratifiedGLTest(RSS)
    res = Any[]
    for i = 1:length(RSS)
        append!(res, 0)
        res[i] = []
        TimerReset()
        for j = 1:length(RSS[i])
            append!(res[i], 0)
            res[i][j] = LScoreCommunity(B, RSS[i][j])[1]
        end
        println(TimerLapValue())
    end
    return res
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