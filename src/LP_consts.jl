using Glob
# Solvers
SOLVER_FN_ADS = 1
SOLVER_LP_ADSS = 2
NUM_SOLVERS = 2
ALL_SOLVERS = [true, true]
LP_SOLVER_ONLY = [false, true]
SOLVER_NAMES = ["FNLA", "LPLAS"]
ERR_TIME_LIMIT = 9999.0 # The time "reported" if exception when solving LP or having no results
TIME_LIMIT = 300.0 # Cutoff time for LP solvers

# Stats to output for LP_algorithm.DoSolveLocalADS (core algorithm for experiment of running ADS problems on FN and LP).
STATS_DS = 1 # result set as a densestSubgraph struct
STATS_OUTPUT_ALPHA = STATS_DS # in the output file "lpcompsets", it is just alpha star value
STATS_EXT_TIME = 2 # external time elapsed
STATS_INT_TIME = 3 # internal time "officially" given
STATS_LNSIZE = 4 # number of explored nodes of the subgraph the last iteration
STATS_LMSIZE = 5 # number of edges of the above subgraph
STATS_ITERS = 6 # number of iterations
STATS_LAST = STATS_ITERS
STATS_OUTPUT_SSIZE = STATS_LAST + 1
STATS_OUTPUT_LAST = STATS_OUTPUT_SSIZE
STATS_NAMES = ["alpha", "ext_time", "int_time", "lnsize", "lmsize", "iters"]
STATS_OUTPUT_NAMES = vcat(STATS_NAMES, ["ssize"])

# I/O of experiment results
RESULT_TYPE_STATS = 1
RESULT_TYPE_SETS = 2
RESULT_TYPE_NAMES = ["lpcompstats", "lpcompsets"]
FOLDER_LP_COMP_RESULTS = "../LPCompResults/"

function GetLPCompResultFileName(dataName::String, solverID::Int, suffixName::String, resultType::Int)
    name = string(dataName, "-", SOLVER_NAMES[solverID == SOLVER_FN_ADS ? SOLVER_FN_ADS : SOLVER_LP_ADSS])
    if length(suffixName) > 0
        name = string(name, "-", suffixName)
    end
    name = string(name, ".", RESULT_TYPE_NAMES[resultType])
    return name
end

function GetLPEvalResultFileName(dataName::String, suffixName::String)
    name = dataName
    if length(suffixName) > 0
        name = string(name, "-", suffixName)
    end
    return name
end

function GetParameterizedLPResultFileNames(dataName::String, solverID::Int64, suffixName::String, resultType::Int)
    return glob(string(FOLDER_LP_COMP_RESULTS, dataName, "-", SOLVER_NAMES[solverID], "-*-*-", suffixName, ".", RESULT_TYPE_NAMES[resultType]), ".")
end
