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
# suffixNames = ["FN100","ADSLsmartL","ADSFsmartL","ADSLSsmartL","ADSFSsmartL"]
# suffixNames = ["FN100","ADSLsmartL","ADSLSsmartL"]
queryID = 11
B = readIN(string(dataName, ".in"))


# Take dataset and query id,
# Generate a csv that contains the union of all result sets.
# Each line represents a node, contains:
# id,orgid,degree,(whether in...)R,ADS,LPA,...,LIAS

# Give R,G,B to specific indexes of resultSets.
# Typically, resultSets contains: ["R","FN100","ADSLsmartL","ADSFsmartL","ADSIsmartL","ADSLSsmartL","ADSFSsmartL","ADSISsmartL"]
# 0 for skipping this color.
indicatorIDs = [3,4,5]  # Means "ADSLsmartL","ADSFsmartL","ADSIsmartL" are for RGB respectively.
# indicatorIDs = [3,4,0]
# indicatorIDs = [3,0,0]

# Result is known
# V: Indicate the single query node, can be -1.
function generateGephiFromAnchor(B::SparseMatrixCSC, dataName::String, V::Int64, suffixNames::Array{String}, queryID::Int64, outputSuffix::String, nodeNames=nothing)
    resultSets = Array{Any}(undef, length(suffixNames) + 1)
    resultSets[1] = readAnchors(dataName, "Baseline")[queryID] # R
    for solverID in eachindex(suffixNames)
        resultSets[solverID + 1] = readCompsets(dataName, solverID, suffixNames[solverID])[queryID]
    end
    unionRes = sort(reduce(union, resultSets))
    df = generateDataframeGephi(B, dataName, V, resultSets, suffixNames, nodeNames)
    CSV.write(string(folderString(FOLDER_CS_LP), join(["cslp", "node", dataName, queryID, outputSuffix], "-"), ".csv"), df, header=true)
    exportGephiEdgelist(B[unionRes, unionRes], string(join(["cslp", "edge", dataName, queryID, outputSuffix], "-"), ".csv"), folderString(FOLDER_CS_LP))
end

function generateDataframeGephi(B::SparseMatrixCSC, dataName::String, V::Int64, resultSets, suffixNames::Array{String}, nodeNames)
    unionRes = sort(reduce(union, resultSets))
    colTypes = vcat(repeat([Int64], 4), repeat([String], 2), repeat([Int64], length(resultSets)))
    colNames = vcat(["id", "orgid", "V", "degree", "Label", "color", "R"], suffixNames)
    df = DataFrame([Vector{t}() for t in colTypes], colNames)
    for id in eachindex(unionRes)
        orgid = unionRes[id]
        isv = orgid == V ? 1 : 0
        degree = GetDegree(B, orgid)
        label = string(degree, isnothing(nodeNames) ? "" : string(", ", nodeNames[orgid]))
        indicator = map(x->orgid in x, resultSets)
        colorInds = Array{Int}(undef, 3)
        for i in 1:3
            colorInds[i] = indicatorIDs[i] > 0 ? indicator[indicatorIDs[i]] : 0
        end
        color = getColorSchemeRGB(colorInds[1], colorInds[2], colorInds[3])
        push!(df, vcat([id, orgid, isv, degree, label, color], indicator))
    end
    return df
end

# A color scheme that takes 4 indicators, the first indicator is considered "black base".
# Each indicator is an array of 0/1s.
function getColorSchemeRGB(r, g, b)
    rc = 0.2 + r * 0.8
    gc = 0.2 + g * 0.8
    bc = 0.2 + b * 0.8
    return rgbToHex(rc, gc, bc)
end
    
    