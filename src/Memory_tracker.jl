# https://stackoverflow.com/questions/38986764/save-variable-name-as-string-in-julia
macro varname(arg)
    string(arg)
end

# This tool manually tracks memory usage of functions.

GLOBAL_memory_dict = Dict{String,Dict{String,Int64}}()
GLOBAL_max_memory_usage = 0
GLOBAL_current_stack = 1 # Depth of current call (among functions of interest only)
GLOBAL_current_stamp = 1 # Number of calls of functions of interest

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

function FunctionKeyName(FunctionName::String, Stamp::Int)
    return string(FunctionName, "-", Stamp)
end

# Call at the start of a function
# stamp = RegisterFunctionStamp()
function RegisterFunctionStamp()
    global GLOBAL_current_stack
    global GLOBAL_current_stamp
    current_stamp = GLOBAL_current_stamp
    GLOBAL_current_stamp += 1
    GLOBAL_current_stack += 1
    return current_stamp
end

function RegisterMemoryItem(FunctionName::String, Stamp::Int, Var::Any, VarName::String)
    global GLOBAL_memory_dict
    functionKey = FunctionKeyName(FunctionName, Stamp)
    if !haskey(GLOBAL_memory_dict, functionKey)
        GLOBAL_memory_dict[functionKey] = Dict{String,Int64}()
    end
    GLOBAL_memory_dict[functionKey][VarName] = Base.summarysize(Var)
    UpdateMaxMemoryUsage()
end

function DeregisterMemoryItem(FunctionName::String, Stamp::Int, VarName::String)
    global GLOBAL_memory_dict
    functionKey = FunctionKeyName(FunctionName, Stamp)
    if haskey(GLOBAL_memory_dict, functionKey)
        delete!(GLOBAL_memory_dict[functionKey], VarName)
    end
end

# Call at the end of the function
function ReclaimFunctionMemoryUsage(FunctionName::String, Stamp::Int)
    global GLOBAL_current_stack
    UpdateMaxMemoryUsage()
    GLOBAL_current_stack -= 1
    delete!(GLOBAL_memory_dict, FunctionKeyName(FunctionName, Stamp))
    # GLOBAL_memory_dict[FunctionKeyName(FunctionName, Stamp)] = Dict{String,Int64}() # If want to leave trace
end

# Take and clear max memory usage value
function PopMaxMemoryUsage(ResetStackAndStamp::Bool=true)
    global GLOBAL_max_memory_usage
    global GLOBAL_current_stack
    global GLOBAL_current_stamp
    memory_usage = GLOBAL_max_memory_usage
    GLOBAL_max_memory_usage = 0
    if ResetStackAndStamp
        GLOBAL_current_stack = 1
        GLOBAL_current_stamp = 1
    end
    return memory_usage
end
