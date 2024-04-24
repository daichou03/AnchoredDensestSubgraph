using HTTP
using JSON
using SparseArrays  # Assuming SparseArrays for demonstration
include("Utils_io.jl")
include("LP_algorithm.jl")

# Placeholder global variable to store the graph
global G = nothing

function ReadGraph(filename::AbstractString, dir::String=DIR_EXAMPLE_SCC)
    return readIN(filename, 1.0, dir)
end

function DoSolveGADS(R, weight_map)
    global G
    return DoSolveLocalADS(SOLVER_LP_ADSS, G, R, false, false, DEFAULT_LP_SOLVER, weight_map).source_nodes
end

function handle_request(request::HTTP.Request)
    try
        data = JSON.parse(String(request.body))
        action = data["action"]
        
        if action == "dummy"
            return HTTP.Response(200, "Dummy action")

        elseif action == "load-graph"
            filename = data["filename"]
            dir = get(data, "dir", DIR_EXAMPLE_SCC)
            global G = ReadGraph(filename, dir)
            return HTTP.Response(200, "Graph loaded successfully from: $(joinpath(dir, filename))")
        
        elseif action == "solve-gads"
            R = Int64.(data["R"])
            weight_map = Float64.(data["weight_map"])
            result = DoSolveGADS(R, weight_map)
            return HTTP.Response(200, JSON.json(result))
        
        else
            return HTTP.Response(400, "Invalid action")
        end
    catch e
        return HTTP.Response(500, "Internal Server Error: $e")
    end
end

HTTP.serve(handle_request, "127.0.0.1", 8080)
