using MAT
include("LP_consts.jl")
include("Utils.jl")
include("Utils_io.jl")
include("CP_FlowSeed.jl")
##########################
# Parameterized GLC test #
##########################
# Similar to Parameterized LP test section in LP_compare_test.jl
TEST_GLC_PENALTY_R = [0.0, 0.01, 0.031, 0.1, 0.31, 1.0, 3.1, 10.0, 31.0, 100.0]
TEST_GLC_EPSILON = [0.01, 0.031, 0.1, 0.31, 1.0, 3.1, 10.0, 31.0, 100.0]


# 20250123: Note that the original function binary searches alpha in the outer loop and S in the inner loop.
function ProcessGLC(B::SparseMatrixCSC, anchors::Array{Array{Int,1},1}, pr::Float64, epsilon::Float64, printInterval=1)
    result_set = Array{Any}(undef, length(anchors))
    prev_time = now()
    for i = 1:length(anchors)
        R = anchors[i]
        result_set[i] = LocalCond(B, R, pr, Int64[], epsilon, true)
        if printInterval > 0 && (now()-prev_time).value / (printInterval * 1000) > 1
            print(string(i, " | ", result_set[i], "\n"))
            prev_time = now()
        end
    end
    return result_set
end


STATS_GLC_COND = 1 # local conductance value
STATS_GLC_EXT_TIME = 2 # external time elapsed
STATS_GLC_LAST = STATS_GLC_EXT_TIME
STATS_GLC_NAMES = ["cond", "ext_time"]
FS_S = 1
FS_COND = 2

function OutputStatsGLC(result_set, dataName::String, suffixName::String="")
    mkpath(FOLDER_LP_COMP_RESULTS)
    io_stats = open(string(FOLDER_LP_COMP_RESULTS, GetLPCompResultFileName(dataName, SOLVER_FN_ADS, suffixName, RESULT_TYPE_STATS)), "w")  # Shouldn't really use SOLVER_FN_ADS though.
    io_sets = open(string(FOLDER_LP_COMP_RESULTS, GetLPCompResultFileName(dataName, SOLVER_FN_ADS, suffixName, RESULT_TYPE_SETS)), "w")
    write(io_stats, string(join(STATS_GLC_NAMES, ","), "\n"))
    for i in 1:length(result_set)
        stats = []
        append!(stats, result_set[i][STATS_GLC_COND][FS_COND])
        for j in (STATS_GLC_COND+1):STATS_GLC_LAST
            append!(stats, result_set[i][j])
        end
        write(io_stats, string(join(map(string, stats),","), "\n"))
        write(io_sets, string(join(map(string, result_set[i][STATS_GLC_COND][FS_S]),","), "\n"))
    end
    close(io_stats)
    close(io_sets)
end


# Similar to LP_compare_test.ProcessAndOutputParameterizedLP
function ProcessAndOutputParameterizedGLC(dataName::String; prRange = TEST_GLC_PENALTY_R, epsRange = TEST_GLC_EPSILON, anchorsType::String="Baseline", suffixName::String="GLCtest", sampleSize::Int=0)
    B = readIN(string(dataName, ".in"))
    anchors = readAnchors(dataName, anchorsType)
    if sampleSize > 0
        anchors = anchors[1:sampleSize]
    end
    for pr in prRange
        for epsilon in epsRange
            println(string("penalty of including R = ", pr, ", epsilon = ", epsilon, ":"))
            result_set = ProcessGLC(B, anchors, pr, epsilon)
            OutputStatsGLC(result_set, dataName, join([pr,epsilon,suffixName], "-"))
        end
    end
end