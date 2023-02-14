using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using Base
using JuMP
include("Helper_io.jl")
include("Graph_utils_yd.jl")
include("Utils.jl")
include("Core_algorithm_yd.jl")

SOLVER_FN_ADS = 1
SOLVER_LP_ADSS = 2
NUM_SOLVERS = 2
ALL_SOLVERS = [true, true]
SOLVER_NAMES = ["FNLA", "LPLAS"]
TIME_LIMIT = 300.0


# Currently support these LP solvers: HiGHS, GLPK, Clp, CDDLib, Gurobi, CPLEX
# Set DEFAULT_LP_SOLVER to change a solver.
# using HiGHS
# DEFAULT_LP_SOLVER = HiGHS

# CPLEX: https://github.com/jump-dev/CPLEX.jl
# Licensed version installed required
# ENV["CPLEX_STUDIO_BINARIES"] = "C:\\Program Files\\IBM\\ILOG\\CPLEX_Studio2211\\cplex\\bin\\x64_win64\\"
# ENV["CPLEX_STUDIO_BINARIES"] = "/opt/ibm/ILOG/CPLEX_Studio201/cplex/bin/x86-64_linux"
# Julia env: "../env/cplex"
using CPLEX
DEFAULT_LP_SOLVER = CPLEX

# Gurobi: https://github.com/jump-dev/Gurobi.jl
# Licensed version installed required
# ENV["GUROBI_HOME"] = "C:\\gurobi1001\\win64\\"
# Julia env: "../env/gurobi"
# using Gurobi
# DEFAULT_LP_SOLVER = Gurobi

# Returns:
# struct:densestSubgraph, time of LP.
function SolveLPDensestSubgraph(B::SparseMatrixCSC, solver=DEFAULT_LP_SOLVER)
    model = SetupLPSolver(DEFAULT_LP_SOLVER)
    edgelist = CSCToEdgeListUndirected(B)
    n = B.n
    m = length(edgelist)
    @variable(model, x[i = 1:n] >= 0)
    @variable(model, y[i = 1:m] >= 0)
    wy = Array{Int}(undef, m)
    @constraint(model, sum(x[i] for i in 1:n) <= 1)
    @objective(model, Max, sum(y))
    for i = 1:m
        u, v = edgelist[i]
        @constraint(model, y[i] <= x[u])
        @constraint(model, y[i] <= x[v])
    end
    optimize!(model)
    return densestSubgraph(objective_value(model), findall(x->value(x)>0, x)), solve_time(model)
end

# No extra output from solver
# Set time limit to TIME_LIMIT
function SetupLPSolver(solver)
    if @isdefined(CDDLib) && solver == CDDLib
        model = Model(solver.Optimizer{Float64})
    else
        model = Model(solver.Optimizer)
    end
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
    end
    return model
end

# Global-LP-ADS#
# Note that "Anchored Densest Subgraph Sharp" means ADS#, which is different from ADS (ADS can't be LP engineered)
function SolveLPAnchoredDensestSubgraphSharp(B::SparseMatrixCSC, R::Vector{Int64}, solver=DEFAULT_LP_SOLVER)
    model = SetupLPSolver(solver)
    edgelist = CSCToEdgeListUndirected(B)
    n = B.n
    m = length(edgelist)
    @variable(model, x[i = 1:n] >= 0)
    @variable(model, y[i = 1:m] >= 0)
    wy = Array{Int}(undef, m)
    @constraint(model, sum(x[i] for i in 1:n) <= 1)
    for i = 1:m
        u, v = edgelist[i]
        if (u in R) || (v in R)
            if (u in R) && (v in R)
                wy[i] = 2
            else
                wy[i] = 1
            end
            @constraint(model, y[i] <= x[u])
            @constraint(model, y[i] <= x[v])
        else
            wy[i] = -1
            @constraint(model, y[i] >= x[u] - x[v])
            @constraint(model, y[i] >= x[v] - x[u])
        end
    end
    @objective(model, Max, sum(map(*, wy, y)))
    optimize!(model)
    return densestSubgraph(objective_value(model), findall(x->value(x)>0, x)), solve_time(model)
end

# Local-LP-ADS#
function SolveLPLocalAnchoredDensestSubgraphSharp(B::SparseMatrixCSC, R::Vector{Int64}, MoreStats::Bool=false, ShowTrace::Bool=false, lpSolver=DEFAULT_LP_SOLVER)
    return DoSolveLocalADS(SOLVER_LP_ADSS, B, R, MoreStats, ShowTrace, lpSolver)
end


function DoSolveLocalADS(Solver::Int, B::SparseMatrixCSC, R::Vector{Int64}, MoreStats::Bool=false, ShowTrace::Bool=false, lpSolver=DEFAULT_LP_SOLVER)
    Expanded = Int64[]
    RSorted = sort(R)
    Frontier = RSorted
    alpha = 0
    S = Int64[]
    SUnion = Int64[]
    L = Int64[]
    total_time = 0
    iters = 0
    while !isempty(Frontier)
        Expanded = union(Expanded, Frontier)
        L = sort(union(L, GetComponentAdjacency(B, Frontier, true))) # GetComponentAdjacency is expensive, doing it incrementally.
        if Solver == SOLVER_FN_ADS
            # TODO: This is unfair, ideally should reconstruct the solution, which contains the time of core algorithms only.
            result_timed = @timed GlobalAnchoredDensestSubgraph(B[L,L], orderedSubsetIndices(L, RSorted))
            result_S, time_taken = result_timed.value, result_timed.time
        elseif Solver == SOLVER_LP_ADSS
            result_S, time_taken = SolveLPAnchoredDensestSubgraphSharp(B[L,L], orderedSubsetIndices(L, RSorted), lpSolver)
        else
            error("Unexpected Solver ID")
        end
        total_time += time_taken
        alpha = result_S.alpha_star
        S = L[result_S.source_nodes]
        if ShowTrace
            println(densestSubgraph(result_S.alpha_star, S))
        end
        SUnion = union(SUnion, S)
        Frontier = setdiff(S, Expanded)
        iters += 1
    end
    
    if MoreStats
        return densestSubgraph(alpha, S), total_time, length(L), iters
    else
        return densestSubgraph(alpha, S)
    end
end
