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
function ProcessAlgorithms(B::SparseMatrixCSC, anchors::Array{Array{Int,1},1}, SolverMask::Vector{Bool}=[true, true], lpWeightMap=DEFAULT_WEIGHT_MAP, printInterval=1)
    statsAlgorithms = []
    for solver_index in 1:length(SolverMask)
        if SolverMask[solver_index]
            result_set = Array{Any}(undef, length(anchors))
            prev_time = now()
            for i = 1:length(anchors)
                R = anchors[i]
                result_set[i] = DoSolveLocalADS(solver_index, B, R, true, false, DEFAULT_LP_SOLVER, lpWeightMap)
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


function ProcessAndOutputAlgorithms(dataName::String, anchorsType="Baseline", SolverMask::Vector{Bool}=ALL_SOLVERS, lpWeightMap=DEFAULT_WEIGHT_MAP, suffixName::String="", sampleSize::Int=0)
    B = readIN(string(dataName, ".in"))
    anchors = readAnchors(dataName, anchorsType)
    if sampleSize > 0
        anchors = anchors[1:sampleSize]
    end
    statsAlgorithms = ProcessAlgorithms(B, anchors, SolverMask, lpWeightMap)
    OutputStatsAlgorithms(statsAlgorithms, dataName, suffixName)
end


dataset_names_1to100m = ["amazon","notredame","digg","citeseer","livemocha","flickr","hyves","yahoo","youtube","google","trec","flixster","dblp","skitter","indian","libimseti","pokec","usaroad","livejournal"]
dataset_names = ["amazon","astroph","brightkite","condmat","dblp","deezer","douban","enron","epinion","fbgov","github","gowalla","grqc","hamster","hepph","hepth","lastfm","livejournal","livemocha","orkut","youtube"]
# 20230228: union of the above two, but excluding libimseti as process is killed when running LP.
dataset_union_to100m = ["amazon", "notredame", "digg", "citeseer", "livemocha", "flickr", "hyves", "yahoo", "youtube", "google", "trec", "flixster", "dblp", "skitter", "indian", "pokec", "usaroad", "livejournal", "astroph", "brightkite", "condmat", "deezer", "douban", "enron", "epinion", "fbgov", "github", "gowalla", "grqc", "hamster", "hepph", "hepth", "lastfm", "orkut"]
dataset_names_lps = ["amazon","notredame","digg","citeseer"] # A
dataset_names_lps = ["livemocha","flickr","hyves","youtube"] # B
dataset_names_lps = ["google","trec","flixster","dblp","skitter"] # C
dataset_names_lps = ["indian","pokec","usaroad","livejournal","orkut"] # E
dataset_names_lps = ["amazon","notredame","digg","citeseer","livemocha","flickr","hyves","youtube","google","trec","flixster","dblp","skitter","indian","pokec","usaroad","livejournal","orkut"]

function BulkProcessAndOutputAlgorithms(dataset_names, SolverMask=ALL_SOLVERS, lpWeightMap=DEFAULT_WEIGHT_MAP, suffixName::String="", sampleSize::Int=0)
    for dataName in dataset_names
        println(string(dataName, ":"))
        proc = @timed ProcessAndOutputAlgorithms(dataName, "Baseline", SolverMask, lpWeightMap, suffixName, sampleSize)
        println(proc.time)
    end
end


############################
# Sensitivity Test: R-Size #
############################

TARGET_SIZES = [8,16,32,64,128,256,512]
function ProcessAndOutputLPFixedSizes(dataName::String, sizes, sampleSize::Int=0)
    B = readIN(string(dataName, ".in"))
    for rsize in sizes
        anchors = readAnchors(dataName, string("fix-", rsize))
        if sampleSize > 0
            anchors = anchors[1:sampleSize]
        end
        statsAlgorithms = ProcessAlgorithms(B, anchors, LP_SOLVER_ONLY)
        OutputStatsAlgorithms(statsAlgorithms, dataName, string("fix-", rsize))
    end
end


#########################
# Parameterized LP test #
#########################

TEST_WAC_RANGE = collect(0:0.2:2)
TEST_WAD_RANGE = [0, -0.0001, -0.00031, -0.001, -0.0031, -0.01, -0.031, -0.1, -0.31, -1]
# w_CC = 1 in parameterized LP tests.
function ProcessAndOutputParameterizedLP(dataName::String; wACRange = TEST_WAC_RANGE, wADRange = TEST_WAD_RANGE, anchorsType::String="Baseline", suffixName::String="", sampleSize::Int=0)
    B = readIN(string(dataName, ".in"))
    anchors = readAnchors(dataName, anchorsType)
    if sampleSize > 0
        anchors = anchors[1:sampleSize]
    end
    for wAC in wACRange
        for wAD in wADRange
            println(string("wAC = ", wAC, ", wAD = ", wAD, ":"))
            lpWeightMap = [2,wAC,0,0,0,0,wAD]
            statsAlgorithms = ProcessAlgorithms(B, anchors, LP_SOLVER_ONLY, lpWeightMap)
            OutputStatsAlgorithms(statsAlgorithms, dataName, join([wAC,abs(wAD),suffixName], "-"))
        end
    end
end


function BulkProcessAndOutputParameterizedLP(dataset_names; wACRange = TEST_WAC_RANGE, wADRange = TEST_WAD_RANGE, anchorsType="Baseline", suffixName::String="", sampleSize::Int=0)
    for dataName in dataset_names
        println(string(dataName, ":"))
        proc = @timed ProcessAndOutputParameterizedLP(dataName; wACRange, wADRange, anchorsType, suffixName, sampleSize)
        println(proc.time)
    end
end


# Halved graphs only ("amazon-H1.in" etc.)
function BulkProcessAndOutputParameterizedLPHalfGraph(dataset_names; wACRange = TEST_WAC_RANGE, wADRange = TEST_WAD_RANGE, anchorsType="Baseline", suffixName::String="", iteration = 5, sampleSize::Int=0)
    halved_dataset_names = String[]
    for dataname in dataset_names
        for i in 1:iteration
            push!(halved_dataset_names, string(dataname, "-H", i))
        end
    end
    BulkProcessAndOutputParameterizedLP(halved_dataset_names; wACRange, wADRange, anchorsType, suffixName, sampleSize)
end

# -------

warmed_up_solver = false

function WarmUpSolvers()
    global warmed_up_solver
    if !warmed_up_solver
        for i= 1:NUM_SOLVERS
            print(string("Warming up solver #", i, "\n"))
            DoSolveLocalADS(i, SAMPLE_GRAPH, SAMPLE_GRAPH_R, false, false, DEFAULT_LP_SOLVER)
        end
        warmed_up_solver = true
        println(string("Done."))
        println(string("DEFAULT_LP_SOLVER = ", DEFAULT_LP_SOLVER))
    end
end

WarmUpSolvers()
