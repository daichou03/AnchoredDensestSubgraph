using CSV
using DataFrames
using StatsBase
include("Utils.jl")
include("LP_consts.jl")

EVAL_DATA_NAME = 1
EVAL_DATA_INDEX = 2
EVAL_ALPHA_EQUAL = 3
EVAL_ALPHA_DIFF = 4
EVAL_EXT_TIME_1 = 5
EVAL_EXT_TIME_2 = 6
EVAL_INT_TIME_1 = 7
EVAL_INT_TIME_2 = 8
EVAL_LN1 = 9
EVAL_LM1 = 10
EVAL_LN2 = 11
EVAL_LM2 = 12
EVAL_SS1 = 13
EVAL_SS2 = 14
EVAL_LAST = EVAL_SS2
EVAL_NAMES = ["data_name", "index", "alpha_equal", "alpha_diff", "ext_time_1", "ext_time_2", "int_time_1", "int_time_2", "ln1", "lm1", "ln2", "lm2", "s_size_1", "s_size_2"]
FOLDER_LP_EVAL_RESULTS = "../LPEvalResults/"

# suffixName: either a string, or an array containing a suffix for each solverID.
function CompareResultSets(dataName::String, suffixName)
    dfs = Array{Any}(undef, 2)
    for solverID in 1:NUM_SOLVERS
        currSuffixName = suffixName isa String ? suffixName : suffixName[solverID]
        dfs[solverID] = DataFrame(CSV.File(string(FOLDER_LP_COMP_RESULTS, GetLPCompResultFileName(dataName, solverID, currSuffixName, RESULT_TYPE_STATS))))
    end
    nrows = min(nrow(dfs[1]), nrow(dfs[2]))
    dataNameColumn = repeat([dataName], nrows)
    dataIndex = collect(1:nrows)
    alphaEquals = Array{Float64}(undef, nrows)
    alphaDiffs = Array{Float64}(undef, nrows)
    timeExt1 = Array{Float64}(undef, nrows)
    timeExt2 = Array{Float64}(undef, nrows)
    timeInt1 = Array{Float64}(undef, nrows)
    timeInt2 = Array{Float64}(undef, nrows)
    lns1 = Array{Float64}(undef, nrows)
    lms1 = Array{Float64}(undef, nrows)
    lns2 = Array{Float64}(undef, nrows)
    lms2 = Array{Float64}(undef, nrows)
    lss1 = Array{Float64}(undef, nrows)
    lss2 = Array{Float64}(undef, nrows)
    for i in 1:nrows
        alphaEquals[i] = almostEqual(dfs[1][i,STATS_NAMES[STATS_OUTPUT_ALPHA]], dfs[2][i,STATS_NAMES[STATS_OUTPUT_ALPHA]]) ? 1 : 0
        alphaDiffs[i] = dfs[2][i,STATS_NAMES[STATS_OUTPUT_ALPHA]] / dfs[1][i,STATS_NAMES[STATS_OUTPUT_ALPHA]]
        timeExt1[i] = dfs[1][i,STATS_NAMES[STATS_EXT_TIME]]
        timeExt2[i] = dfs[2][i,STATS_NAMES[STATS_EXT_TIME]]
        timeInt1[i] = dfs[1][i,STATS_NAMES[STATS_INT_TIME]]
        timeInt2[i] = dfs[2][i,STATS_NAMES[STATS_INT_TIME]]
        lns1[i] = dfs[1][i,STATS_NAMES[STATS_LNSIZE]]
        lms1[i] = dfs[1][i,STATS_NAMES[STATS_LMSIZE]]
        lns2[i] = dfs[2][i,STATS_NAMES[STATS_LNSIZE]]
        lms2[i] = dfs[2][i,STATS_NAMES[STATS_LMSIZE]]
        lss1[i] = dfs[1][i,STATS_NAMES[STATS_OUTPUT_SSIZE]]
        lss2[i] = dfs[2][i,STATS_NAMES[STATS_OUTPUT_SSIZE]]
    end
    return dataNameColumn, dataIndex, alphaEquals, alphaDiffs, timeExt1, timeExt2, timeInt1, timeInt2, lns1, lms1, lns2, lms2, lss1, lss2
end

function OutputCompareResults(evalResults, dataName::String, suffixName::String)
    mkpath(FOLDER_LP_EVAL_RESULTS)
    io = open(string(FOLDER_LP_EVAL_RESULTS, GetLPEvalResultFileName(dataName, suffixName)), "w")
    write(io, string(join(EVAL_NAMES, ","), "\n"))
    for i in 1:length(evalResults[1])
        stats = Array{Any}(undef, EVAL_LAST)
        for j in 1:EVAL_LAST
            stats[j] = evalResults[j][i]
        end
        write(io, string(join(stats, ","), "\n"))
    end
    close(io)
end

function CompareAndOutputResultSets(dataName::String, suffixName)
    evalResults = CompareResultSets(dataName, suffixName)
    OutputCompareResults(evalResults, dataName, suffixName isa String ? suffixName : suffixName[1])
end

function BulkCompareAndOutputResultSets(dataNames, suffixName)
    for dataName in dataNames
        CompareAndOutputResultSets(dataName, suffixName)
    end
end

function BulkCompareAndOutputConcatenatedResultSets(dataNames, concatName, suffixName)
    resultCon = Array{Any}(undef, EVAL_LAST)
    for i in 1:EVAL_LAST
        resultCon[i] = []
    end
    for dataName in dataNames
        result = CompareResultSets(dataName, suffixName)
        for i in 1:EVAL_LAST
            append!(resultCon[i], result[i])
        end
    end
    OutputCompareResults(resultCon, concatName, suffixName)
end


# TODO:
# Also output LM, LN (have minimal difference, take both first), regression on these.
# Optimizaitons