using CSV
using DataFrames
using StatsBase
include("Utils.jl")
include("Utils_io.jl")
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

FILENAME_EMPTY_LPCOMPSTATS = "empty.lpcompstats"
# 20230323
# Take ADS (Flow Network), ADSL, ADSF together for example.
# For example, assume .lpcompstats has name like "amazon-FNLA-FN100.lpcompstats", "google-LPLAS-ADSL100C.lpcompstats", etc.
# Sample dataNames:
# ["amazon", "google"]
# Sample suffixNames - one for each model to compare. Note assuming the first name is for FNLA, and the second name onwards are all for LPLAS.
# suffixNames = ["FN100", "ADSL100C", "ADSF100C"]
# If a file is not found, assume the task couldn't be completed and fall back to use FILENAME_EMPTY_LPCOMPSTATS as data.
function CompareMultipleModelResultSets(dataName::String, suffixNames::Array{String}, getRatio::Bool=false)
    df1 = Nothing
    means = Array{Any}(undef, length(suffixNames))
    for solverID in 1:length(suffixNames)
        filename = GetLPCompResultFileName(dataName, solverID == 1 ? SOLVER_FN_ADS : SOLVER_LP_ADSS, suffixNames[solverID], RESULT_TYPE_STATS)
        fileFound = isfile(string(FOLDER_LP_COMP_RESULTS, filename))
        if fileFound
            df = DataFrame(CSV.File(string(FOLDER_LP_COMP_RESULTS, filename)))
        else
            df = DataFrame(CSV.File(string(FOLDER_LP_COMP_RESULTS, FILENAME_EMPTY_LPCOMPSTATS)))
        end
        if getRatio
            if solverID == 1
                df1 = copy(df)
            end
            if fileFound
                df = df ./ df1
            end
        end
        means[solverID] = [mean(df[!, column]) for column in names(df)]
    end
    return means
end

# Example: given dataMeans like:
# Amazon
# alpha,ext_time,int_time
# 1,0.011,0.005 <- ADS
# 2,0.024,0.011 <- ADSL
# 3,0.024,0.011 <- ADSF
# DBLP
# alpha,ext_time,int_time
# 4,0.011,0.005 <- ADS
# 5,0.024,0.011 <- ADSL
# 6,0.024,0.011 <- ADSF
# Output file with name average-suffix-alpha (also for ext_time and int_time):
# data,ADS,ADSL,ADSF
# amazon,1,2,3
# DBLP,4,5,6

suffixNames = ["FN100","ADSL100C","ADSF100C","ADSI100C","ADSLS100C","ADSFS100C","ADSIS100C"]
function OutputMultipleModelResultSets(dataNames::Array{String}, suffixNames::Array{String}, outputSuffix::String, getRatio::Bool=false)
    dataMeans = Array{Any}(undef, length(dataNames))
    for dataID in eachindex(dataNames)
        dataMeans[dataID] = CompareMultipleModelResultSets(dataNames[dataID], suffixNames, getRatio)
    end
    columnNames = names(DataFrame(CSV.File(string(FOLDER_LP_COMP_RESULTS, GetLPCompResultFileName(
        dataNames[1], SOLVER_FN_ADS, suffixNames[1], RESULT_TYPE_STATS)))))
    col_names = vcat(["dataName"], suffixNames)
    col_types = vcat(String, repeat([Float64], length(suffixNames)))

    mkpath(FOLDER_LP_EVAL_RESULTS)
    for columnID in eachindex(columnNames)
        df = DataFrame([Vector{t}() for t in col_types], col_names)
        for dataID in 1:length(dataNames)
            push!(df, vcat(dataNames[dataID], [row[columnID] for row in dataMeans[dataID]]))
        end
        CSV.write(string(folderString(FOLDER_LP_EVAL_RESULTS), join([getRatio ? "ratio" : "average", outputSuffix, columnNames[columnID]], "-")), df, header=true)
    end
end


# For F1 score
function readCompsets(DataName::AbstractString, SolverID, SuffixName::String, EmptyIfNotFound::Bool=false, SubDirName::String="")
    filename = string(folderString(FOLDER_LP_COMP_RESULTS, SubDirName), GetLPCompResultFileName(DataName, SolverID, SuffixName, RESULT_TYPE_SETS)) 
    if EmptyIfNotFound && !isfile(filename)
        return []
    end
    results = []
    for rawline in eachline(filename)
        line = strip(rawline)
        if length(line) > 0
            line = split(line, ",")
            line = map(x->parse(Int, x), line)
            append!(results, [line])
        end
    end
    return results
end


function CompareMultipleModelF1score(dataName::String, suffixNames::Array{String})
    results = []
    means = Array{Any}(undef, length(suffixNames))
    anchors = readAnchors(dataName, "Baseline")
    for algID in 1:length(suffixNames)
        solverID = algID == 1 ? SOLVER_FN_ADS : SOLVER_LP_ADSS
        results = readCompsets(dataName, solverID, suffixNames[algID], true)
        resultSetLength = min(length(anchors), length(results))
        means[solverID] = mean(map(i->f1score(anchors[i], results[i]), 1:resultSetLength))
    end
    return means
end

function OutputMultipleModelF1score(dataNames::Array{String}, suffixNames::Array{String}, outputSuffix::String)
    dataMeans = Array{Any}(undef, length(dataNames))
    for dataID in eachindex(dataNames)
        dataMeans[dataID] = CompareMultipleModelF1score(dataNames[dataID], suffixNames)
    end
    col_names = vcat(["dataName"], suffixNames)
    col_types = vcat(String, repeat([Float64], length(suffixNames)))

    mkpath(FOLDER_LP_EVAL_RESULTS)
    df = DataFrame([Vector{t}() for t in col_types], col_names)
    for dataID in 1:length(dataNames)
        push!(df, vcat(dataNames[dataID], [row[columnID] for row in dataMeans[dataID]]))
    end
    CSV.write(string(folderString(FOLDER_LP_EVAL_RESULTS), join(["average", outputSuffix, "f1score"], "-")), df, header=true)
end
