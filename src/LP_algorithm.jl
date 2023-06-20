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


DEFAULT_WEIGHT_MAP = WEIGHT_MAP_ADSL


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
OPT_SMARTL = true # If true, the starting working graph does not include edges between neighbours of R. This reduces starting |E(L)| at the cost of more iterations. Experiment shows it is beneficial to LP but not so for FN.
# Strongly-local version that integrates both flow network and linear programming solution.
function DoSolveLocalADS(Solver::Int, B::SparseMatrixCSC, R::Vector{Int64}, MoreStats::Bool=false, ShowTrace::Bool=false, lpSolver=DEFAULT_LP_SOLVER)
    R = sort(R)
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
    # Expanded = Int64[]
    # Frontier = R
    alpha = 0.0
    S = Int64[] # Result set for each iteration
    C = sort(OPT_SMARTL ? R : GetComponentAdjacency(B, R, true)) # Inner nodes - all edges
    P = R # Result nodes new to SRUnion. Start with non-empty.
    SRUnion = sort(union(S,R))
    F = sort(setdiff(GetComponentAdjacency(B, P, false), C)) # Frontier nodes - nodes are in the working graph, but edges between two F nodes do not.
    L = [B[C,C] B[C,F]; B[F,C] spzeros(length(F),length(F))] # Working subgraph
    int_time = 0
    ext_time = 0
    iters = 0
    inducedDS = GlobalDensestSubgraph(B[R,R])
    if inducedDS.alpha_star >= 1 # If not, return empty results early. This is necessary to have both maximal and minimal algorithm behave similarly in this edge case.
        while !isempty(P)
            overdensedMask = map(v->(GetDegree(B,v)>=GetVolume(B,R)), [C;F])
            if Solver == SOLVER_FN_ADS
                result_timed = @timed ImprovedGlobalAnchoredDensestSubgraphSetFlow(L, orderedSubsetIndices([C;F], R), overdensedMask, inducedDS)
                result_S, ext_time_taken = result_timed.value, result_timed.time
                # Take ext_time as int_time for now.
                int_time_taken = ext_time_taken
            elseif Solver == SOLVER_LP_ADSS
                mip_set = [findfirst([C;F] .== v) for v in (length(S) > 0 ? S : R[inducedDS.source_nodes])]
                if WeightIsADSIX(DEFAULT_WEIGHT_MAP) # IGA optimization is not correct for ADSIX.
                    overdensedMask = nothing
                end
                result_timed = @timed SolveLPAnchoredDensestSubgraphGeneric(L, orderedSubsetIndices([C;F], R), DEFAULT_WEIGHT_MAP, overdensedMask, mip_set, lpSolver)
                ext_time_taken = result_timed.time
                result_S, int_time_taken = result_timed.value
            else
                error("Unexpected Solver ID")
            end
            alpha = result_S.alpha_star
            int_time += int_time_taken
            ext_time += ext_time_taken
            S = [C;F][result_S.source_nodes]
            C = sort(union(C, S))
            P = sort(setdiff(S, SRUnion))
            SRUnion = sort(union(S,R))
            F = sort(setdiff(union(F, GetComponentAdjacency(B, P, false)), C))
            L = [B[C,C] B[C,F]; B[F,C] spzeros(length(F),length(F))]
            iters += 1
            if ShowTrace
                println(join([densestSubgraph(result_S.alpha_star, S), ext_time, int_time, L.n, nnz(L)÷2, iters, length(S), length(C), length(P), length(SRUnion), length(F)], " | "))
            end
            # ADSFX weights always one-shot
            if WeightIsADSFX(DEFAULT_WEIGHT_MAP)
                break
            end
        end
    end
    
    if MoreStats
        # See LP_consts.STATS_NAMES
        return densestSubgraph(alpha, S), ext_time, int_time, L.n, nnz(L)÷2, iters
    else
        return densestSubgraph(alpha, S)
    end
end

