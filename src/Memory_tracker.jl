# https://stackoverflow.com/questions/38986764/save-variable-name-as-string-in-julia
macro varname(arg)
    string(arg)
end

GLOBAL_memory_dict = Dict{String,Dict{String,Int64}}()
GLOBAL_max_memory_usage = 0

function UpdateMaxMemoryUsage()
    global GLOBAL_memory_dict
    global GLOBAL_max_memory_usage
    sum_memory_usage = 0
    for (func_name, func_memory) in GLOBAL_memory_dict
        for (var_name, var_memory) in func_memory
            sum_memory_usage += var_memory
        end
    end
    GLOBAL_max_memory_usage = max(GLOBAL_max_memory_usage, sum_memory_usage)
end

# Call after each major assignment
function RegisterMemoryItem(FunctionName::String, Var::Any, VarName::String)
    global GLOBAL_memory_dict
    if !haskey(GLOBAL_memory_dict, FunctionName)
        GLOBAL_memory_dict[FunctionName] = Dict{String,Int64}()
    end
    GLOBAL_memory_dict[FunctionName][VarName] = Base.summarysize(Var)
    UpdateMaxMemoryUsage()
end

function DeregisterMemoryItem(FunctionName::String, VarName::String)
    global GLOBAL_memory_dict
    if haskey(GLOBAL_memory_dict, FunctionName) && haskey(GLOBAL_memory_dict[FunctionName], VarName)
        delete!(GLOBAL_memory_dict[FunctionName], VarName)
    end
end

# Call before return from a function
function ReclaimFunctionMemoryUsage(FunctionName::String)
    UpdateMaxMemoryUsage()
    GLOBAL_memory_dict[FunctionName] = Dict{String,Int64}()
end

# Take and clear max memory usage value
function PopMaxMemoryUsage()
    global GLOBAL_max_memory_usage
    memory_usage = GLOBAL_max_memory_usage
    GLOBAL_max_memory_usage = 0
    return memory_usage
end

# function BulkRegisterMemoryItems(FunctionName::String, Vars::Array{Any, 1}, VarNames::Array{String, 1})
#     for i in 1:length(Vars)
#         RegisterMemoryItem(FunctionName, Vars[i], VarNames[i])
#     end
# end
