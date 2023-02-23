using CSV
using DataFrames
using StatsBase
include("Utils.jl")
include("LP_consts.jl")

########################################
# Compare LA and Local-LP-ADS# results #
########################################

# Simply count number of result sets that are equal.
# function CompareResultSets(dataName::String, suffixName::String="")
#     hit, miss = 0, 0
#     dfs = Array{Any}(undef, 2)
#     for solverID in 1:NUM_SOLVERS
#         dfs[solverID] = DataFrame(CSV.File(string(FOLDER_LP_COMP_RESULTS, GetLPCompResultFileName(dataName, solverID, suffixName, RESULT_TYPE_STATS))))
#     end
#     for i in 1:length(dfs[1].alpha)
#         almostEqual(dfs[1].alpha[i], dfs[2].alpha[i]) ? hit += 1 : miss += 1
#     end
#     alphaDiff = mean(dfs[2].alpha) / mean(dfs[1].alpha)
#     time1, time2 = mean(dfs[1].time), mean(dfs[2].time)
#     return join(map(string, [hit, miss, alphaDiff, time1, time2]),",")
# end


# suffixName: either a string, or an array containing a suffix for each solverID.
function CompareResultSets(dataName::String, suffixName)
    hit, miss = 0, 0
    dfs = Array{Any}(undef, 2)
    for solverID in 1:NUM_SOLVERS
        currSuffixName = suffixName isa String ? suffixName : [suffixName[solverID]]
        dfs[solverID] = DataFrame(CSV.File(string(FOLDER_LP_COMP_RESULTS, GetLPCompResultFileName(dataName, solverID, currSuffixName, RESULT_TYPE_STATS))))
    end
    nrows = min(nrow(dfs[1]), nrow(dfs[2]))
    alphaEquals = Array{Float64}(undef, nrows)
    alphaDiffs = Array{Float64}(undef, nrows)
    timeExtDiffs = Array{Float64}(undef, nrows)
    timeIntDiffs = Array{Float64}(undef, nrows)
    for i in 1:length(dfs[1][:,STATS_NAMES[STATS_ALPHA]])nrow(dfs[1])
        alphaEquals[i] = almostEqual(dfs[1][i,STATS_NAMES[STATS_ALPHA]], dfs[2][i,STATS_NAMES[STATS_ALPHA]]) ? 1 : 0
        alphaDiffs[i] = dfs[2][i,STATS_NAMES[STATS_ALPHA]] / dfs[1][i,STATS_NAMES[STATS_ALPHA]]
        timeExtDiffs[i] = dfs[2][i,STATS_NAMES[STATS_EXT_TIME]] / dfs[1][i,STATS_NAMES[STATS_EXT_TIME]]
        timeIntDiffs[i] = dfs[2][i,STATS_NAMES[STATS_INT_TIME]] / dfs[1][i,STATS_NAMES[STATS_INT_TIME]]
    end
    return alphaEquals, alphaDiffs, timeExtDiffs, timeIntDiffs
end