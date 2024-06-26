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
#include("CP_GreedyL.jl")
include("CS_Amazon.jl")

# 20211023: Old integration of startified amazon only tests. May not work without fixing.
# Amazon Stratified tests

function BulkReportCommunity(B::SparseMatrixCSC, Rs::Any, Ss::Any, TestName::String, AlgName::String)
    folder = string(CS_AMAZON_FOLDER, "Report/", TestName, "/EV-", AlgName, "/")
    mkpath(folder)
    for i in 1:length(Rs)
        io = open(string(folder,i,".txt"), "w")
        for j in 1:length(Rs[i])
            write(io, string(ReportCommunity(B, Rs[i][j], Ss[i][j]), "\n"))
        end
        close(io)
    end
    BulkReportRSize(Rs, TestName)
end

function BulkReportRSize(Rs::Any, TestName::String)
    folder = string(CS_AMAZON_FOLDER, "Report/", TestName, "/R/")
    mkpath(folder)
    for i in 1:length(Rs)
        io = open(string(folder,i,".txt"), "w")
        for j in 1:length(Rs[i])
            write(io, string(length(Rs[i][j]), "\n"))
        end
        close(io)
    end
end

####################
# Integrate report #
####################

NUM_REPORTS = 41

function IntegrateCSReport(TestName::String, NumReports::Int64=NUM_REPORTS)
    statsAlgs = ReadCSReport(TestName, NumReports)
    rLengths = ReadCSRLength(TestName, NumReports)
    output_folder = string(CS_AMAZON_FOLDER, "ReportIntegrated/", TestName, "/")
    mkpath(output_folder)
    OutputCSMetrics(statsAlgs, rLengths, output_folder, NumReports)

    calcStatsAlgss = []
    stat_f1score = ReadCSF1score(TestName, NumReports)
    OutputCSCalculatedMetric(stat_f1score, CALCULATED_METRICS[IND_F1SCORE], output_folder, NumReports)
    append!(calcStatsAlgss, [stat_f1score])

    OutputCSAggregated(statsAlgs, calcStatsAlgss, output_folder, NumReports)
end

function ReadCSReport(TestName::String, NumReports::Int64=NUM_REPORTS)
    statsAlgs = []
    for i_alg in 1:length(ALG_REPORT_NAMES)
        folder = string(CS_AMAZON_FOLDER, "Report/", TestName, "/", ALG_REPORT_NAMES[i_alg], "/")
        statsDegs = []
        for i_deg in 1:NumReports
            io = open(string(folder,i_deg,".txt"))
            stats = zeros(length(REPORT_METRICS))
            count = 0
            while !eof(io)
                stat = ImputeNaNs(map(x->parse(Float64, x), split(readline(io), "|")))
                count += 1
                for j in 1:length(REPORT_METRICS)
                    stats[j] += stat[j]
                end
            end
            close(io)
            append!(statsDegs, 0)
            statsDegs[i_deg] = map(x->x/count, stats)
        end
        append!(statsAlgs, 0)
        statsAlgs[i_alg] = statsDegs
    end
    return statsAlgs
end

function ReadCSRLength(TestName::String, NumReports::Int64=NUM_REPORTS)
    folder = string(CS_AMAZON_FOLDER, "Report/", TestName, "/R/")
    rLengths = []
    for i_deg in 1:NumReports
        io = open(string(folder,i_deg,".txt"))
        rLength = 0
        count = 0
        while !eof(io)
            rLength += parse(Int64, readline(io))
            count += 1
        end
        close(io)
        append!(rLengths, 0)
        rLengths[i_deg] = map(x->x/count, rLength)
    end
    return rLengths
end

# For each metrics output data
function OutputCSMetrics(statsAlgs, rLengths, output_folder, NumReports::Int64=NUM_REPORTS)
    for i_metric in 1:length(REPORT_METRICS)
        io = open(string(output_folder, REPORT_METRICS[i_metric], ".txt"), "w")
        for i_deg in 1:NumReports
            line = []
            for i_alg in 1:length(ALG_REPORT_NAMES)
                append!(line, 0)
                line[i_alg] = statsAlgs[i_alg][i_deg][i_metric]
            end
            # Append length of R
            if i_metric == 1
                append!(line, 0)
                line[4] = rLengths[i_deg]
            end
            write(io, string(join(line, " "), "\n"))
        end
        close(io)
    end
end

function ReadCSF1score(TestName::String, NumReports::Int64=NUM_REPORTS)
    calcStatAlgs = []
    for i_alg in 1:length(ALG_REPORT_NAMES)
        folder = string(CS_AMAZON_FOLDER, "Report/", TestName, "/", ALG_REPORT_NAMES[i_alg], "/")
        statsDegs = []
        for i_deg in 1:NumReports
            io = open(string(folder,i_deg,".txt"))
            f1 = 0.0
            count = 0
            while !eof(io)
                stat = ImputeNaNs(map(x->parse(Float64, x), split(readline(io), "|")))
                f1 += f1score(stat[IND_RINS], stat[IND_SINR])
                count += 1
            end
            close(io)
            append!(statsDegs, 0)
            statsDegs[i_deg] = f1 / count
        end
        append!(calcStatAlgs, 0)
        calcStatAlgs[i_alg] = statsDegs
    end
    return calcStatAlgs
end

function OutputCSCalculatedMetric(calcStatAlgs, calcStatName, output_folder, NumReports::Int64=NUM_REPORTS)
    io = open(string(output_folder, calcStatName, ".txt"), "w")
    for i_deg in 1:NumReports
        line = []
        for i_alg in 1:length(ALG_REPORT_NAMES)
            append!(line, 0)
            line[i_alg] = calcStatAlgs[i_alg][i_deg]
        end
        write(io, string(join(line, " "), "\n"))
    end
    close(io)
end

function OutputCSF1score(statsAlgs, output_folder, NumReports::Int64=NUM_REPORTS)
    io = open(string(output_folder, "f1score.txt"), "w")
    for i_deg in 1:NumReports
        line = []
        for i_alg in 1:length(ALG_REPORT_NAMES)
            append!(line, 0)
            line[i_alg] = statsAlgs[i_alg][i_deg]
        end
        write(io, string(join(line, " "), "\n"))
    end
    close(io)
end

function OutputCSAggregated(statsAlgs, calcStatsAlgss, output_folder, NumReports::Int64=NUM_REPORTS)
    io = open(string(output_folder, "aggregated.txt"), "w")
    for i_metric in 1:length(REPORT_METRICS)
        line = [REPORT_METRICS[i_metric]]
        for i_alg in 1:length(ALG_REPORT_NAMES)
            stat = 0.0
            for i_deg in 1:NumReports
                stat += statsAlgs[i_alg][i_deg][i_metric] / NumReports
            end
            append!(line, [string(stat)])
        end
        write(io, string(join(line, " "), "\n"))
    end
    for i_metric in 1:length(calcStatsAlgss)
        line = [CALCULATED_METRICS[i_metric]]
        for i_alg in 1:length(ALG_REPORT_NAMES)
            stat = 0.0
            for i_deg in 1:NumReports
                stat += calcStatsAlgss[i_metric][i_alg][i_deg] / NumReports
            end
            append!(line, [string(stat)])
        end
        write(io, string(join(line, " "), "\n"))
    end
    close(io)
end