using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using StatsBase
using Random
using Base
using DataFrames
using CSV
using Dates
include("maxflow.jl")
include("Helper_io.jl")
include("Graph_utils_yd.jl")
include("Core_algorithm_yd.jl")
include("LP_algorithm.jl")
include("Test_utils_yd.jl")
include("Utils.jl")


STATS_DS = 1
STATS_TIME = 2
STATS_LSIZE = 3
STATS_ITERS = 4
STATS_LAST = STATS_ITERS
STATS_NAMES = ["alpha", "time", "lsize", "iters"]

RESULT_TYPE_STATS = 1
RESULT_TYPE_SETS = 2
RESULT_TYPE_NAMES = ["lpcompstats", "lpcompsets"]

FOLDER_LP_COMP_RESULTS = "../LPCompResults/"

############################
# Run LA and Local-LP-ADS# #
############################
function ProcessAlgorithms(B::SparseMatrixCSC, anchors::Array{Array{Int,1},1}, SolverMask::Vector{Bool}=[true, true], printInterval=-1)
    statsAlgorithms = []
    for solver_index in 1:length(SolverMask)
        if SolverMask[solver_index]
            result_set = Array{Any}(undef, length(anchors))
            prev_time = now()
            for i = 1:length(anchors)
                R = anchors[i]
                result_set[i] = DoSolveLocalADS(solver_index, B, R, true, false, DEFAULT_LP_SOLVER)
                if printInterval > 0 && (now()-prev_time).value / (printInterval * 1000) > 1
                    print(string(i, " | ", result_set[i], "\n"))
                    prev_time = now()
                end
            end
            append!(statsAlgorithms, [result_set])
        else
            append!(statsAlgorithms, [[]])
        end
    end
    return statsAlgorithms
end


function OutputStatsAlgorithms(statsAlgorithms, dataName::String, suffixName::String="")
    for solver_id = 1:NUM_SOLVERS
        if length(statsAlgorithms[solver_id]) > 0
            mkpath(FOLDER_LP_COMP_RESULTS)
            io_stats = open(string(FOLDER_LP_COMP_RESULTS, GetLPCompResultFileName(dataName, solver_id, suffixName, RESULT_TYPE_STATS)), "w")
            io_sets = open(string(FOLDER_LP_COMP_RESULTS, GetLPCompResultFileName(dataName, solver_id, suffixName, RESULT_TYPE_SETS)), "w")
            write(io_stats, string(join(STATS_NAMES, ","), "\n"))
            for i in 1:length(statsAlgorithms[solver_id])
                stats = []
                statsAlgorithms[solver_id][i]
                append!(stats, statsAlgorithms[solver_id][i][STATS_DS].alpha_star)
                for j in STATS_TIME:STATS_LAST
                    append!(stats, statsAlgorithms[solver_id][i][j])
                end
                write(io_stats, string(join(map(string, stats),","), "\n"))
                write(io_sets, string(join(map(string, statsAlgorithms[solver_id][i][STATS_DS].source_nodes),","), "\n"))
            end
            close(io_stats)
            close(io_sets)
        end
    end
end


function ProcessAndOutputAlgorithms(dataName::String, SolverMask::Vector{Bool}=ALL_SOLVERS, suffixName::String="", sampleSize::Int=0)
    B = readIN(string(dataName, ".in"))
    anchors = readAnchors(dataName, "Baseline")
    if sampleSize > 0
        anchors = anchors[1:sampleSize]
    end
    statsAlgorithms = ProcessAlgorithms(B, anchors, SolverMask)
    OutputStatsAlgorithms(statsAlgorithms, dataName, suffixName)
end


dataset_names = ["amazon","astroph","brightkite","condmat","dblp","deezer","douban","enron","epinion","fbgov","github","gowalla","grqc","hamster","hepph","hepth","lastfm","livejournal","livemocha","orkut","youtube"]

function BulkProcessAndOutputAlgorithms(dataset_names, suffixName::String="", sampleSize::Int=0)
    for dataName in dataset_names
        print(string(dataName, ":"))
        proc = @timed ProcessAndOutputAlgorithms(dataName, ALL_SOLVERS, suffixName, sampleSize)
        print(proc.time)
    end
end


function GetLPCompResultFileName(dataName::String, solverID::Int, suffixName::String, resultType::Int)
    name = string(dataName, "-", SOLVER_NAMES[solverID])
    if length(suffixName) > 0
        name = string(name, "-", suffixName)
    end
    name = string(name, ".", RESULT_TYPE_NAMES[resultType])
    return name
end

########################################
# Compare LA and Local-LP-ADS# results #
########################################

# Simply count number of result sets that are equal.
function CompareResultSets(dataName::String, suffixName::String="")
    hit, miss = 0, 0
    dfs = Array{Any}(undef, 2)
    for solverID in 1:NUM_SOLVERS
        dfs[solverID] = DataFrame(CSV.File(string(FOLDER_LP_COMP_RESULTS, GetLPCompResultFileName(dataName, solverID, suffixName, RESULT_TYPE_STATS))))
    end
    for i in 1:length(dfs[1].alpha)
        almostEqual(dfs[1].alpha[i], dfs[2].alpha[i]) ? hit += 1 : miss += 1
    end
    alphaDiff = mean(dfs[2].alpha) / mean(dfs[1].alpha)
    time1, time2 = mean(dfs[1].time), mean(dfs[2].time)
    return join(map(string, [hit, miss, alphaDiff, time1, time2]),",")
end


warmed_up_solver = false

function WarmUpSolvers()
    global warmed_up_solver
    if !warmed_up_solver
        print("Warming up solvers...")
        for i= 1:NUM_SOLVERS
            DoSolveLocalADS(i, SAMPLE_GRAPH, SAMPLE_GRAPH_R, DEFAULT_LP_SOLVER)
        end
        warmed_up_solver = true
        print("Done.")
    end
end

WarmUpSolvers()

# using GLPK
# DEFAULT_LP_SOLVER = GLPK
# dataName = "astroph"
# ProcessAndOutputAlgorithms(dataName, [false, true], "GLPK100", 100)