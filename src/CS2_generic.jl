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
suffixNames = ["FN100","ADSLsmartL","ADSFsmartL","ADSLSsmartL","ADSFSsmartL"]
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
# indicatorIDs = [3,0,4]

# Result is known
function generateGephiFromAnchor(B::SparseMatrixCSC, dataName::String, suffixNames::Array{String}, queryID::Int64, outputSuffix::String)
    resultSets = Array{Any}(undef, length(suffixNames) + 1)
    resultSets[1] = readAnchors(dataName, "Baseline")[queryID] # R
    for solverID in eachindex(suffixNames)
        resultSets[solverID + 1] = readCompsets(dataName, solverID, suffixNames[solverID])[queryID]
    end
    return generateGephi(B, dataName, -1, resultSets, suffixNames, outputSuffix)
end

function generateGephi(B::SparseMatrixCSC, dataName::String, V::Int64, resultSets, suffixNames::Array{String}, outputSuffix::String)
    unionRes = sort(reduce(union, resultSets))
    colTypes = vcat(repeat([Int64], 4), repeat([String], 2), repeat([Int64], length(resultSets)))
    colNames = vcat(["id", "orgid", "V", "degree", "Label", "color", "R"], suffixNames)
    df = DataFrame([Vector{t}() for t in colTypes], colNames)
    for id in eachindex(unionRes)
        orgid = unionRes[id]
        isv = orgid == V ? 1 : 0
        degree = GetDegree(B, orgid)
        label = string(degree)
        indicator = map(x->orgid in x, resultSets)
        colorInds = Array{Int}(undef, 3)
        for i in 1:3
            colorInds[i] = indicatorIDs[i] > 0 ? indicator[indicatorIDs[i]] : 0
        end
        color = getColorSchemeRGB(colorInds[1], colorInds[2], colorInds[3])
        push!(df, vcat([id, orgid, isv, degree, label, color], indicator))
    end
    CSV.write(string(folderString(FOLDER_CS_LP), join(["cslp", "node", dataName, queryID, outputSuffix], "-"), ".csv"), df, header=true)
    exportGephiEdgelist(B[unionRes, unionRes], string(join(["cslp", "edge", dataName, queryID, outputSuffix], "-"), ".csv"), folderString(FOLDER_CS_LP))
end

# A color scheme that takes 4 indicators, the first indicator is considered "black base".
# Each indicator is an array of 0/1s.
function getColorSchemeRGB(r, g, b)
    rc = 0.5 + r * 0.5
    gc = 0.5 + g * 0.5
    bc = 0.5 + b * 0.5
    return rgbToHex(rc, gc, bc)
end
    
    