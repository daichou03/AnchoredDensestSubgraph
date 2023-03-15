using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using Random
using Base
include("Utils_io.jl")
include("Utils_graph.jl")
include("Utils_warmup.jl")
include("Utils.jl")
include("maxflow.jl")
include("Core_algorithm_yd.jl")
include("LP_consts.jl")
include("LP_algorithm.jl")


############################
# Run LA and Local-LP-ADS# #
############################
function ProcessAlgorithms(B::SparseMatrixCSC, anchors::Array{Array{Int,1},1}, SolverMask::Vector{Bool}=[true, true], printInterval=1)
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
            write(io_stats, string(join(STATS_OUTPUT_NAMES, ","), "\n"))
            for i in 1:length(statsAlgorithms[solver_id])
                stats = []
                statsAlgorithms[solver_id][i]
                append!(stats, statsAlgorithms[solver_id][i][STATS_DS].alpha_star)
                for j in (STATS_DS+1):STATS_LAST
                    append!(stats, statsAlgorithms[solver_id][i][j])
                end
                append!(stats, length(statsAlgorithms[solver_id][i][STATS_DS].source_nodes))
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

dataset_names_1to100m = ["amazon","notredame","digg","citeseer","livemocha","flickr","hyves","yahoo","youtube","google","trec","flixster","dblp","skitter","indian","libimseti","pokec","usaroad","livejournal"]
dataset_names = ["amazon","astroph","brightkite","condmat","dblp","deezer","douban","enron","epinion","fbgov","github","gowalla","grqc","hamster","hepph","hepth","lastfm","livejournal","livemocha","orkut","youtube"]
# 20230228: union of the above two, but excluding libimseti as process is killed when running LP.
dataset_union_to100m = ["amazon", "notredame", "digg", "citeseer", "livemocha", "flickr", "hyves", "yahoo", "youtube", "google", "trec", "flixster", "dblp", "skitter", "indian", "pokec", "usaroad", "livejournal", "astroph", "brightkite", "condmat", "deezer", "douban", "enron", "epinion", "fbgov", "github", "gowalla", "grqc", "hamster", "hepph", "hepth", "lastfm", "orkut"]

function BulkProcessAndOutputAlgorithms(dataset_names, SolverMask=ALL_SOLVERS, suffixName::String="", sampleSize::Int=0)
    for dataName in dataset_names
        println(string(dataName, ":"))
        proc = @timed ProcessAndOutputAlgorithms(dataName, SolverMask, suffixName, sampleSize)
        println(proc.time)
    end
end


warmed_up_solver = false

function WarmUpSolvers()
    global warmed_up_solver
    if !warmed_up_solver
        for i= 1:NUM_SOLVERS
            print(string("Warming up solver #", i, "\n"))
            DoSolveLocalADS(i, SAMPLE_GRAPH, SAMPLE_GRAPH_R, false, false, DEFAULT_LP_SOLVER)
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
# dataName = "livejournal"
# ProcessAndOutputAlgorithms(dataName, [false, true], "GLPK100", 100)