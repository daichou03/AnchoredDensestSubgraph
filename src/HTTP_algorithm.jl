using HTTP
using JSON
using SparseArrays  # Assuming SparseArrays for demonstration
include("Utils_io.jl")
include("LP_algorithm.jl")

# Placeholder global variable to store the graph
global G = nothing
global sync_graph_name = ""

function LoadGraph(filename::AbstractString, dir::String=DIR_EXAMPLE_SCC)
    global G, sync_graph_name
    G = readIN(filename, 1.0, dir)
    sync_graph_name = ""
end

function LoadGraphSync(dataName::AbstractString)
    global G, sync_graph_name
    if sync_graph_name != dataName
        G = readIN("csc.in", 1.0, folderString(DIR_EXAMPLE_DEMO_SYNC, dataName))
        sync_graph_name = dataName
    end
end


function DoSolveGADS(R, weight_map)
    global G
    return DoSolveLocalADS(SOLVER_LP_ADSS, G, R, false, false, DEFAULT_LP_SOLVER, weight_map).source_nodes
end

function handle_request(request::HTTP.Request)
    try
        data = JSON.parse(String(request.body))
        action = data["action"]
        
        println(action)
        
        if action == "dummy"
            return HTTP.Response(200, "Dummy action")

        elseif action == "load-graph"
            filename = data["filename"]
            dir = get(data, "dir", DIR_EXAMPLE_SCC)
            LoadGraph(filename, dir)
            return HTTP.Response(200, "Graph loaded successfully from: $(joinpath(dir, filename))")

        elseif action == "load-graph-sync"
            dataname = data["dataname"]
            LoadGraphSync(dataname)
            return HTTP.Response(200, "Synchronized data $(dataname) loaded successfully")
        
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
