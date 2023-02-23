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
include("CS_Simple.jl")
include("Collect_results.jl")

FS_PENALTY_R = 0.0
FS_EPSILON = 1.0

function GetDensity(B::SparseMatrixCSC, S::Vector{Int64})
    return GetVolume(B[S,S]) / length(S)
end

function GetAnchoredDensity(B::SparseMatrixCSC, R::Vector{Int64}, S::Vector{Int64})
    S_in_R = intersect(S, R)
    S_out_R = setdiff(S, S_in_R)
    return (GetVolume(B[S,S]) - GetVolume(B, S_out_R)) / length(S)
end

function GetConductance(B::SparseMatrixCSC, S::Vector{Int64})
    return 1 - GetVolume(B[S,S]) / GetVolume(B, S)
end

function GetLocalConductance(B::SparseMatrixCSC, R::Vector{Int64}, S::Vector{Int64})
    O_R = GetVolume(B, intersect(R, S)) - FS_EPSILON * GetVolume(B, setdiff(S, R)) - FS_PENALTY_R * GetVolume(B, setdiff(R, S))
    return O_R > 0 ? (GetVolume(B, S) - GetVolume(B[S,S])) / O_R : Inf
end

function GetLScore(B::SparseMatrixCSC, S::Vector{Int64})
    return LScore(B, S)
end

function GetPropRinS(R::Vector{Int64}, S::Vector{Int64})
    return 1 - length(setdiff(R, S)) / length(R)
end

# function GetPropSoutR(R::Vector{Int64}, S::Vector{Int64})
#     return length(setdiff(S, R)) / length(S)
# end

function GetPropSinR(R::Vector{Int64}, S::Vector{Int64})
    return 1 - length(setdiff(S, R)) / length(S)
end

function ReportCommunity(B::SparseMatrixCSC, R::Vector{Int64}, S::Vector{Int64})
    return join([length(S),
        GetDensity(B,S),
        GetAnchoredDensity(B,R,S),
        GetConductance(B,S),
        GetLocalConductance(B,R,S),
        #GetLScore(B,S),
        GetPropRinS(R,S),
        GetPropSinR(R,S)], "|")
end

###############
# Bulk Report #
###############

function ReportCommunitySimple(B::SparseMatrixCSC, Rs, Ss, Times, Spaces, DataName::String, AlgName::String)
    folder = folderString(CS_SIMPLE_FOLDER, DataName, "Report")
    mkpath(folder)
    io = open(string(folder, AlgName, ".txt"), "w")
    for j in 1:length(Rs)
        write(io, string(join([Times[j], Spaces[j] / 1000000, ReportCommunity(B, Rs[j], Ss[j])], "|"), "\n"))
    end
    close(io)
end

function BulkReportCommunitySimple()
    for dataName in SIMPLE_TEST_DATA_NAMES
        println(string("Now reporting ", dataName, ":"))
        B = readIN(string(dataName, ".in"))
        vs, rs = ImportSimpleRs(dataName)
        for algName in ALGORITHM_NAMES
            ss, times, spaces = ImportSimpleResults(dataName, algName)
            # Can choose to truncate MRW's result set here
            if algName == "MRW"
                # ss = map(i->ss[i][1:min(length(ss[i]), length(rs[i]))], 1:length(ss)) # Size of R
                ss = map(s->s[1:min(length(s), 15)], ss) # Size = 15
            end
            ReportCommunitySimple(B, rs, ss, times, spaces, dataName, algName)
        end
    end
end

# Fill up spaces. Temporary, later output space too.
function TinkSpaceIntoReport(Spaces, DataName::String, AlgName::String)
    folder = folderString(CS_SIMPLE_FOLDER, DataName, "Report")
    mkpath(folder)
    io = open(string(folder, AlgName, ".txt"))
    arrs = []
    for j in 1:length(Spaces)
        arr = split(readline(io), "|")
        space = string(Spaces[j] / 1000000)
        if length(arr) == 8
            arr = [arr[1];space;arr[2:end]]
        else
            arr[2] = space
        end
        arr = join(arr, "|")
        push!(arrs, arr)
    end
    close(io)
    io = open(string(folder, AlgName, ".txt"), "w")
    for j in 1:length(arrs)
        write(io, string(arrs[j], "\n"))
    end
    close(io)
end

function TinkTimeIntoReport(Times, DataName::String, AlgName::String)
    folder = folderString(CS_SIMPLE_FOLDER, DataName, "Report")
    mkpath(folder)
    io = open(string(folder, AlgName, ".txt"))
    arrs = []
    for j in 1:length(Times)
        arr = split(readline(io), "|")
        time = string(Times[j])
        if length(arr) == 8
            arr = [time;arr[2:end]]
        else
            arr[1] = time
        end
        arr = join(arr, "|")
        push!(arrs, arr)
    end
    close(io)
    io = open(string(folder, AlgName, ".txt"), "w")
    for j in 1:length(arrs)
        write(io, string(arrs[j], "\n"))
    end
    close(io)
end

function BulkTinkSpaceIntoReport()
    for dataName in SIMPLE_TEST_DATA_NAMES
        for algName in ALGORITHM_NAMES
            ss, times, spaces = ImportSimpleResults(dataName, algName)
            TinkSpaceIntoReport(spaces, dataName, algName)
        end
    end
end


####################
# Integrate report #
####################

ALG_REPORT_NAMES = ["EV-LA", "EV-MRW", "EV-FS"]
REPORT_METRICS = ["time", "space", "length", "density", "rsdensity", "conductance", "lconductance", "rins", "sinr"]
IND_LENGTH = 3
IND_RINS = 8
IND_SINR = 9
CALCULATED_METRICS = ["f1score"]
IND_F1SCORE = 1

IMPUTE_VALUES = [1000.0, 0.0, 0.0, 0.0, 1.0, 999999.0, 0.0, 0.0]

function ImputeNaNs(Values::Vector{Float64}, ImputeValues = IMPUTE_VALUES)
    ret = copy(Values)
    for i in 1:length(ret)
        if isnan(ret[i]) || isinf(ret[i])
            ret[i] = ImputeValues[i]
        end
    end
    return ret
end

function IntegrateSimpleCSReport()
    dataStats, calcDataStats = ReadSimpleCSReport()
    rLengths = ReadSimpleCSRLength()
    output_folder = folderString(PERFORMANCE_REPORTS_INTEGRATED_DIR, "cs_simple")
    mkpath(output_folder)
    OutputSimpleCSMetrics(output_folder, dataStats, rLengths)
    OutputSimpleCSCalculatedMetrics(output_folder, calcDataStats)
end

# Data -> Algorithm -> Average Stats
function ReadSimpleCSReport()
    dataStats = []
    calcDataStats = []
    for dataName in SIMPLE_TEST_DATA_NAMES
        algStats = []
        calcAlgStats = []
        for algName in ALGORITHM_NAMES
            io = open(string(folderString(CS_SIMPLE_FOLDER, dataName, "Report"), algName, ".txt"))
            stats = zeros(length(REPORT_METRICS))
            calcStats = zeros(length(CALCULATED_METRICS))
            count = 0
            while !eof(io)
                stat = ImputeNaNs(map(x->parse(Float64, x), split(readline(io), "|")))
                count += 1
                for j in 1:length(REPORT_METRICS)
                    stats[j] += stat[j]
                end
                for j in 1:length(CALCULATED_METRICS)
                    calcStats[j] += CalculatedMetrics(stat, j)
                end
            end
            close(io)
            push!(algStats, map(x->x/count, stats))
            push!(calcAlgStats, map(x->x/count, calcStats))
        end
        push!(dataStats, algStats)
        push!(calcDataStats, calcAlgStats)
    end
    return dataStats, calcDataStats
end

function ReadSimpleCSRLength()
    rLengths = []
    for dataName in SIMPLE_TEST_DATA_NAMES
        io = open(string(folderString(CS_SIMPLE_FOLDER, dataName, "R"), "r.txt"))
        rLength = 0
        count = 0
        while !eof(io)
            rLength += length(split(readline(io), ","))
            count += 1
        end
        close(io)
        push!(rLengths, rLength / count)
    end
    return rLengths
end

# For each metrics output data
function OutputSimpleCSMetrics(folder, dataStats, rLengths)
    for i_metric in 1:length(REPORT_METRICS)
        io = open(string(folder, REPORT_METRICS[i_metric], ".txt"), "w")
        for i_data in 1:length(SIMPLE_TEST_DATA_NAMES)
            line = []
            push!(line, SIMPLE_TEST_DATA_NAMES[i_data])
            for i_alg in 1:length(ALG_REPORT_NAMES)
                push!(line, dataStats[i_data][i_alg][i_metric])
            end
            # Append length of R for result length
            if i_metric == IND_LENGTH
                push!(line, rLengths[i_data])
            end
            write(io, string(join(line, " "), "\n"))
        end
        close(io)
    end
end

function OutputSimpleCSCalculatedMetrics(folder, calcDataStats)
    for i_metric in 1:length(CALCULATED_METRICS)
        io = open(string(folder, CALCULATED_METRICS[i_metric], ".txt"), "w")
        for i_data in 1:length(SIMPLE_TEST_DATA_NAMES)
            line = []
            push!(line, SIMPLE_TEST_DATA_NAMES[i_data])
            for i_alg in 1:length(ALG_REPORT_NAMES)
                push!(line, calcDataStats[i_data][i_alg][i_metric])
            end
            write(io, string(join(line, " "), "\n"))
        end
        close(io)
    end
end

function f1score(p, r)
    if p == 0 && r == 0
        return 0
    end
    return 2 * p * r / (p + r)
end

function CalculatedMetrics(stat, index)
    if index == IND_F1SCORE
        return f1score(stat[IND_RINS], stat[IND_SINR])
    end
    throw(string("Unexpected index of calculated metric index: ", index))
end
