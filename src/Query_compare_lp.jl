using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using StatsBase # TODO: To install
using Random
using Base
include("maxflow.jl")
include("Helper_io.jl")
include("Graph_utils_yd.jl")
include("Core_algorithm_yd.jl")
include("LP_algorithm.jl")
include("Test_utils_yd.jl")
include("Utils.jl")


# Comparing LA and Local-LP-ADS#
function ProcessAlgorithms(B::SparseMatrixCSC, anchors::Array{Array{Int,1},1}, SolverMask::Vector{Bool}=[true, true])
    statsAlgorithms = []
    for solver_index in 1:length(SolverMask)
        if SolverMask[solver_index]
            result_set = Array{Any}(undef, length(anchors))
            for i = 1:length(anchors)
                R = anchors[i]
                result_set[i] = DoSolveLocalADS(solver_index, B, R, true, false)
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
            mkpath("../LPCompResults")
            filename = string("../LPCompResults/", dataName, "-", SOLVER_NAMES[solver_id])
            if length(suffixName) > 0
                filename = string(filename, "-", suffixName)
            end
            io_stats = open(string(filename, ".lpcompstats"), "w")
            io_sets = open(string(filename, ".lpcompsets"), "w")
            write(io_stats, string(join(STATS_NAMES,","), "\n"))
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


function ProcessAndOutputAlgorithms(dataName::String, SolverMask::Vector{Bool}=[true, true], suffixName::String="")
    B = readIN(string(dataName, ".in"))
    anchors = readAnchors(dataName, "Baseline")
    statsAlgorithms = ProcessAlgorithms(B, anchors, SolverMask)
    OutputStatsAlgorithms(statsAlgorithms, dataName, suffixName)
end


dataset_names = ["amazon","astroph","brightkite","condmat","dblp","deezer","douban","enron","epinion","fbgov","github","gowalla","grqc","hamster","hepph","hepth","lastfm","livejournal","livemocha","orkut","youtube"]

function BulkProcessAndOutputAlgorithms(dataset_names, suffixName::String="")
    for dataName in dataset_names
        proc = @timed ProcessAndOutputAlgorithms(dataName, ALL_SOLVERS, suffixName)
        print(string(proc, ": ", proc.time))
    end
end


warmed_up_solver = false

function WarmUpSolvers()
    global warmed_up_solver
    if !warmed_up_solver
        print("Warming up solvers...")
        for i= 1:NUM_SOLVERS
            DoSolveLocalADS(i,SAMPLE_GRAPH,SAMPLE_GRAPH_R)
        end
        warmed_up_solver = true
        print("Done.")
    end
end

WarmUpSolvers()
