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

    colTypes = vcat(repeat([Int64], 3), repeat([String], 2), repeat([Int64], length(results)))
    colNames = vcat(["id", "orgid", "degree", "Label", "color", "R"], suffixNames)
    df = DataFrame([Vector{t}() for t in colTypes], colNames)
    for id in eachindex(unionRes)
        orgid = unionRes[id]
        degree = GetDegree(B, orgid)
        label = string(degree)
        indicator = map(x->orgid in x, results)
        color = getColorSchemeKRGB(indicator[1], indicator[3], indicator[4], indicator[5])
        push!(df, vcat([id, orgid, degree, label, color], indicator))
    end
    CSV.write(string(folderString(FOLDER_CS_LP), join(["cslp", "node", dataName, queryID], "-"), ".csv"), df, header=true)
    exportGephiEdgelist(B[unionRes, unionRes], string(join(["cslp", "edge", dataName, queryID], "-"), ".csv"), folderString(FOLDER_CS_LP))
end

# A color scheme that takes 4 indicators, the first indicator is considered "black base".
# Each indicator is an array of 0/1s.
function getColorSchemeKRGB(k, r, g, b)
    rc = 0.5 + r * 0.5
    gc = 0.5 + g * 0.5
    bc = 0.5 + b * 0.5
    return rgbToHex(rc, gc, bc)
end
    
    