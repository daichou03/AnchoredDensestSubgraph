using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using Base
using JuMP
include("Utils_io.jl")
include("Utils_graph.jl")
include("Utils.jl")
include("Core_algorithm_yd.jl")
include("LP_consts.jl")
include("LP_load_solver.jl")


# Returns:
# struct:densestSubgraph, time of LP.
# TODO: Need to prove if calling local version of this could lead to a sub-optimal local result.
function SolveLPDensestSubgraph(B::SparseMatrixCSC, solver=DEFAULT_LP_SOLVER)
    model = SetupLPSolver(solver)
    edgelist = CSCToEdgeListUndirected(B)
    n = B.n
    m = length(edgelist)
    @variable(model, x[i = 1:n] >= 0)
    @variable(model, y[i = 1:m] >= 0)
    @constraint(model, sum(x[i] for i in 1:n) <= 1)
    @objective(model, Max, 2*sum(y))
    for i = 1:m
        u, v = edgelist[i]
        @constraint(model, y[i] <= x[u])
        @constraint(model, y[i] <= x[v])
    end
    try
        optimize!(model)
        return densestSubgraph(objective_value(model), findall(x->value(x)>0, x)), solve_time(model)
    catch y
        println("Exception: ", y)
        return EMPTY_DENSEST_SUBGRAPH, ERR_TIME_LIMIT
    end
end

PARAM_THREADS = -1 # -1 means do not force set number of threads
# No extra output from solver
# Set time limit to TIME_LIMIT
function SetupLPSolver(solver)
    # Init
    if @isdefined(CDDLib) && solver == CDDLib
        model = Model(solver.Optimizer{Float64})
    else
        model = Model(solver.Optimizer)
    end
    # Slient, time limit
    if @isdefined(HiGHS) && solver == HiGHS
        set_optimizer_attribute(model, "log_to_console", false)
        set_optimizer_attribute(model, "time_limit", TIME_LIMIT)
    elseif @isdefined(GLPK) && solver == GLPK
        set_optimizer_attribute(model, "msg_lev", GLPK.GLP_MSG_OFF)
        set_optimizer_attribute(model, "tm_lim", TIME_LIMIT * 1000)
    elseif @isdefined(Clp) && solver == Clp
        set_optimizer_attribute(model, "LogLevel", 0)
        set_optimizer_attribute(model, "MaximumSeconds", TIME_LIMIT)
    elseif @isdefined(Gurobi) && solver == Gurobi
        set_optimizer_attribute(model, "TimeLimit", TIME_LIMIT)
        set_optimizer_attribute(model, "LogToConsole", 0)
    else
        set_silent(model)
        set_time_limit_sec(model, TIME_LIMIT)
    end
    # Number of threads (CPLEX, HiGHS only)
    if PARAM_THREADS >= 0
        if @isdefined(CPLEX) && solver == CPLEX
            set_optimizer_attribute(model, "CPX_PARAM_THREADS", PARAM_THREADS)
        elseif @isdefined(HiGHS) && solver == HiGHS
            set_optimizer_attribute(model, "threads", PARAM_THREADS)
        end
    end
    return model
end

OPT_IGA = true # IGA optimization: Set overdensed node's value to 0
OPT_MIP = false # MIP start optimization
OPT_FEASIBILITY = false # Currently CPLEX only: Higher feasibility to match problem
OPT_DUAL = false # Currently CPLEX only: Solve dual problem instead
OPT_NOPRESOLVE = false # Currently CPLEX only: Do not presolve
OPT_BARRIER = false # Currently CPLEX only: Use barrier algorithm


# Weight maps. 7 positions stand for ... accordingly:
WEIGHT_IND_RSXRS = 1    # R∩S x R∩S
WEIGHT_IND_RSXS = 2     # R∩S x S
WEIGHT_IND_RSXR = 3     # R∩S x R
WEIGHT_IND_RSXE = 4     # R∩S x ∅
WEIGHT_IND_SXS = 5      # S x S
WEIGHT_IND_SXR = 6      # S x R
WEIGHT_IND_SXE = 7      # S x ∅

WEIGHT_MAP_DS = [2,2,0,0,2,0,0] # (Global) Densest Subgraph. Runnable on LP, but not strongly local nor optimal on local algorithm.
# WEIGHT_MAP_ADS = [2,1,0,0,0,-1,-1] # Anchored Densest Subgraph for flow network. Note that LP algorithm won't work for this weight map.
WEIGHT_MAP_ADSL = [2,1,0,0,0,0,-1] # Anchored Densest Subgraph Linear (previously called ADS Sharp, the first LP variation researched). 20230324: proven to be strongly local nor optimal on local algorithm.
WEIGHT_MAP_ADSF = [2,1,0,0,0,0,0] # Anchored Densest Subgraph Fast. 20230324: proven to be not strongly local nor optimal on local algorithm.
WEIGHT_MAP_ADSI = [2,1,0,0,1,0,0] # Anchored Densest Subgraph Intense. 20230324: proven to be NOT strongly local nor optimal on local algorithm.
WEIGHT_MAP_ADSLS = [2,2,0,0,0,0,-1] # ADSL, but weight of R∩S x S is 2. 20230324: All conclusions on locality for ADSXS follows the non-S version.
WEIGHT_MAP_ADSFS = [2,2,0,0,0,0,0] # ADSF, but weight of R∩S x S is 2
WEIGHT_MAP_ADSIS = [2,2,0,0,1,0,0] # ADSI, but weight of R∩S x S is 2

function WeightIsADSFX(weightMap)
    return all(x -> x > 0, weightMap[1:2]) && all(x -> almostEqual(x, 0), weightMap[3:7])
end

DEFAULT_WEIGHT_MAP = WEIGHT_MAP_ADSL


# EIR0_TYPE_X: X = |e∩R|.
EIR0_TYPE_ZERO = 0
EIR0_TYPE_POSITIVE = 1
EIR0_TYPE_NEGATIVE = -1

WEIGHT_FEATURE_EIR = [2,1,2,1,0,1,0] # Labels EIR type for each edge type, not an actual WEIGHT_MAP.

# Checks whether the weight values for |e∩R| = 0 edges are all zero, all non-negative or all non-positive.
# Note this is only a problem for |e∩R| = 0 edges, as |e∩R| = 1 or 2 edges must be all non-negative.
function FindWeightMapEIR0Type(WeightMap)
    weightMapEIR0 = findall(WEIGHT_FEATURE_EIR .== 0)
    if all(x -> x == 0, WeightMap[weightMapEIR0])
        return EIR0_TYPE_ZERO
    elseif all(x -> x >= 0, WeightMap[weightMapEIR0])
        return EIR0_TYPE_POSITIVE # All non-negative actually
    elseif all(x -> x <= 0, WeightMap[weightMapEIR0])
        return EIR0_TYPE_NEGATIVE # All non-positive actually
    else
        throw(ArgumentError("All values in WeightMap under the same EIR type must be either non-positive or non-negative"))
    end
end


# Note that you can only assume this algorithm to work as expected with ONLY weight maps defined above,  
# as they all meet some *SPECIFIC* assumptions (see paper).
function SolveLPAnchoredDensestSubgraphGeneric(B::SparseMatrixCSC, R::Vector{Int64}, WeightMap = WEIGHT_MAP_ADSL, OverdensedMask=nothing, MIP=nothing, solver=DEFAULT_LP_SOLVER)
    model = SetupLPSolver(solver)
    edgelist = CSCToEdgeListUndirected(B)
    n = B.n
    m = length(edgelist)
    @variable(model, x[i = 1:n] >= 0)
    @variable(model, y[i = 1:m] >= 0)
    wy = Array{Float64}(undef, m)
    @constraint(model, sum(x[i] for i in 1:n) <= 1)

    # Optimizations
    if OPT_IGA && !isnothing(OverdensedMask)
        for i = 1:n
            if OverdensedMask[i]
                @constraint(model, x[i] == 0)
            end
        end
    end
    
    if OPT_MIP && !isnothing(MIP)
        for i = 1:n
            if i in MIP
                set_start_value(x[i], 1 / length(MIP))
            else
                set_start_value(x[i], 0)
            end
        end            
    end

    if OPT_FEASIBILITY
        set_optimizer_attribute(model, "CPX_PARAM_EPRHS", 0.005)
    end

    if OPT_DUAL
        set_optimizer_attribute(model, "CPXPARAM_Simplex_Pricing", CPX_ALG_DUAL)
    end

    if OPT_NOPRESOLVE
        set_optimizer_attribute(model, "CPX_PARAM_PREIND", 0)
    end

    if OPT_BARRIER
        set_optimizer_attribute(model, "CPX_PARAM_LPMETHOD", 4)
    end

    # Weights and objective
    eir0type = FindWeightMapEIR0Type(WeightMap)
    for i = 1:m
        u, v = edgelist[i]
        eir = (u in R ? 1 : 0) + (v in R ? 1 : 0) # |e∩R|
        if eir == 0 && eir0type == EIR0_TYPE_NEGATIVE
            wy[i] = WeightMap[WEIGHT_IND_SXE]
            @constraint(model, y[i] >= x[u] - x[v])
            @constraint(model, y[i] >= x[v] - x[u])
        elseif eir == 0 && eir0type == EIR0_TYPE_ZERO
            wy[i] = 0
            @constraint(model, y[i] == 0)
        else
            if eir == 2
                wy[i] = WeightMap[WEIGHT_IND_RSXRS]
            elseif eir == 1
                wy[i] = WeightMap[WEIGHT_IND_RSXS]
            else # eir == 0 && eir0type == EIR0_TYPE_POSITIVE
                wy[i] = WeightMap[WEIGHT_IND_SXS]
            end
            @constraint(model, y[i] <= x[u])
            @constraint(model, y[i] <= x[v])
        end
    end
    @objective(model, Max, sum(map(*, wy, y)))
    # 20230221: Possible reason of exceptions:
    # Something's wrong when optimize!()
    # Get no result
    try
        optimize!(model)
        return densestSubgraph(objective_value(model), findall(x->value(x)>0, x)), solve_time(model)
    catch y
        println("Exception: ", y)
        return EMPTY_DENSEST_SUBGRAPH, ERR_TIME_LIMIT
    end
end


OPT_FORCE_GLOBAL = false # Force running global version rather than running strongly-local
# Strongly-local version that integrates both flow network and linear programming solution.
function DoSolveLocalADS(Solver::Int, B::SparseMatrixCSC, R::Vector{Int64}, MoreStats::Bool=false, ShowTrace::Bool=false, lpSolver=DEFAULT_LP_SOLVER)
    if OPT_FORCE_GLOBAL && Solver == SOLVER_LP_ADSS
        mip_set = GlobalDensestSubgraph(B[R,R]).source_nodes 
        result_timed = @timed SolveLPAnchoredDensestSubgraphGeneric(B, R, DEFAULT_WEIGHT_MAP, nothing, mip_set, lpSolver)
        ext_time_taken = result_timed.time
        result_S, int_time_taken = result_timed.value
        alpha = result_S.alpha_star
        if MoreStats
            # See LP_consts.STATS_NAMES
            return result_S, ext_time_taken, int_time_taken, B.n, nnz(B)÷2, 1
        else
            return result_S
        end
    end
    Expanded = Int64[]
    RSorted = sort(R)
    Frontier = RSorted
    alpha = 0
    S = Int64[]
    SUnion = Int64[]
    L = Int64[]
    int_time = 0
    ext_time = 0
    iters = 0
    inducedDS = GlobalDensestSubgraph(B[R,R])
    while !isempty(Frontier)
        Expanded = union(Expanded, Frontier)
        L = sort(union(L, GetComponentAdjacency(B, Frontier, true))) # GetComponentAdjacency is expensive, doing it incrementally.
        overdensedMask = map(v->(GetDegree(B,v)>=GetVolume(B,R)), L)
        if Solver == SOLVER_FN_ADS
            result_timed = @timed ImprovedGlobalAnchoredDensestSubgraphSetFlow(B[L,L], orderedSubsetIndices(L, RSorted), overdensedMask, inducedDS)
            result_S, ext_time_taken = result_timed.value, result_timed.time
            # Take ext_time as int_time for now.
            int_time_taken = ext_time_taken
        elseif Solver == SOLVER_LP_ADSS
            mip_set = length(S) > 0 ? [findfirst(L .== v) for v in S] : inducedDS.source_nodes
            result_timed = @timed SolveLPAnchoredDensestSubgraphGeneric(B[L,L], orderedSubsetIndices(L, RSorted), DEFAULT_WEIGHT_MAP, overdensedMask, mip_set, lpSolver)
            ext_time_taken = result_timed.time
            result_S, int_time_taken = result_timed.value
        else
            error("Unexpected Solver ID")
        end
        int_time += int_time_taken
        ext_time += ext_time_taken
        alpha = result_S.alpha_star
        S = L[result_S.source_nodes]
        if ShowTrace
            println(densestSubgraph(result_S.alpha_star, S))
        end
        SUnion = union(SUnion, S)
        Frontier = setdiff(S, Expanded)
        iters += 1
        # ADSFX weights always one-shot
        if WeightIsADSFX(DEFAULT_WEIGHT_MAP)
            break
        end
    end
    
    if MoreStats
        # See LP_consts.STATS_NAMES
        return densestSubgraph(alpha, S), ext_time, int_time, length(L), nnz(B[L,L])÷2, iters
    else
        return densestSubgraph(alpha, S)
    end
end

