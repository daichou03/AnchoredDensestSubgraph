using CSV
using DataFrames
using StatsBase
include("../utils/Utils.jl")
include("../GADS/lp_consts.jl")
include("../GADS/lp_evaluation.jl")
include("generic_v2.jl")
include("dblp.jl")

DATA_NAME = "csdblp"
suffixNames = ["FNcsdblp","ADSLcsdblp","ADSFcsdblp","ADSIcsdblp","ADSLScsdblp","ADSFScsdblp","ADSIScsdblp"]
# suffixNames = ["FNcsdblp","ADSLcsdblp","ADSFcsdblp","ADSLScsdblp","ADSFScsdblp"]
indicatorIDs = [3,4,5]

# suffixNames = ["FNcsdblp","ADSLcsdblp","ADSFcsdblp","ADSLScsdblp","ADSFScsdblp"]
# indicatorIDs = [3,4,0]

# suffixNames = ["FNcsdblp","ADSLcsdblp","ADSLScsdblp"]
# indicatorIDs = [3,0,0]
FOLDER_CS_DBLP_CANDIDATE_LP = folderString(FOLDER_CS_DBLP, "candidate", "LP")
Vs = map(x->parse(Int,x), readlines(string(FOLDER_CS_DBLP_CANDIDATE_LP, "Vs.txt")))

function generateGephiFromAnchorCSDBLP(queryID::Int64, suffixNames::Array{String}, outputSuffix::String)
    # B, allNames from CS_DBLP.jl
    return generateGephiFromAnchor(B, DATA_NAME, Vs[queryID], suffixNames, queryID, outputSuffix, allNames)
end
