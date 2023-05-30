using CSV
using DataFrames
using StatsBase
include("Utils.jl")
include("Utils_io.jl")
include("LP_consts.jl")
include("Utils_graph.jl")
include("LP_evaluation.jl")

# CS2 is for ADS revisited - LP_....
FOLDER_CS_LP = "../CaseStudy/LP/"


dataName = "flixster"
suffixNames = ["FN100","ADSLsmartL","ADSFsmartL","ADSIsmartL","ADSLSsmartL","ADSFSsmartL","ADSISsmartL"]
queryID = 11
B = readIN(string(dataName, ".in"))


# Take dataset and query id,
# Generate a csv that contains the union of all result sets.
# Each line represents a node, contains:
# id,orgid,degree,(whether in...)R,ADS,LPA,...,LIAS

function generateGephi(B::SparseMatrixCSC, dataName::String, suffixNames::Array{String}, queryID::Int64)
    results = Array{Any}(undef, length(suffixNames) + 1)
    results[1] = readAnchors(dataName, "Baseline")[queryID] # R
    for solverID in eachindex(suffixNames)
        results[solverID + 1] = readCompsets(dataName, solverID, suffixNames[solverID])[queryID]
    end
    unionRes = sort(reduce(union, results))

    colTypes = repeat([Int64], length(suffixNames) + 4)
    colNames = vcat(["id", "orgid", "degree", "R"], suffixNames)
    df = DataFrame([Vector{t}() for t in colTypes], colNames)
    for id in eachindex(unionRes)
        orgid = unionRes[id]
        push!(df, vcat([id, orgid, GetDegree(B, orgid)], map(x->orgid in x, results)))
    end
    CSV.write(string(folderString(FOLDER_CS_LP), join(["cslp", "node", dataName, queryID], "-"), ".csv"), df, header=true)
    exportGephiEdgelist(B[unionRes, unionRes], string(join(["cslp", "edge", dataName, queryID], "-"), ".csv"), folderString(FOLDER_CS_LP))
end
