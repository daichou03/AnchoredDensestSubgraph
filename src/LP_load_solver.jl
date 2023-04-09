using JuMP
include("Utils_io.jl")
include("Utils_graph.jl")
include("Utils.jl")
include("LP_consts.jl")

# Currently support these LP solvers: HiGHS, GLPK, Clp, CDDLib, Gurobi, CPLEX
# Set DEFAULT_LP_SOLVER to change a solver.
# if !@isdefined(DEFAULT_LP_SOLVER)
#     using HiGHS
#     DEFAULT_LP_SOLVER = HiGHS
# end

# CPLEX: https://github.com/jump-dev/CPLEX.jl
# Licensed version installed required
# ENV["CPLEX_STUDIO_BINARIES"] = "C:\\Program Files\\IBM\\ILOG\\CPLEX_Studio2211\\cplex\\bin\\x64_win64\\"
# ENV["CPLEX_STUDIO_BINARIES"] = "/opt/ibm/ILOG/CPLEX_Studio201/cplex/bin/x86-64_linux"
# Julia env: "../env/cplex"
if !@isdefined(DEFAULT_LP_SOLVER)
    using CPLEX
    DEFAULT_LP_SOLVER = CPLEX
end

# Gurobi: https://github.com/jump-dev/Gurobi.jl
# Licensed version installed required
# ENV["GUROBI_HOME"] = "C:\\gurobi1001\\win64\\"
# Julia env: "../env/gurobi"
# if !@isdefined(DEFAULT_LP_SOLVER)
#     using Gurobi
#     DEFAULT_LP_SOLVER = Gurobi
# end