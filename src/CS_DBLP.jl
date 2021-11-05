using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using StatsBase
using Random
using Base
include("maxflow.jl")
include("Helper_io.jl")
include("Graph_utils_yd.jl")
include("Core_algorithm_yd.jl")
include("Test_utils_yd.jl")
include("Utils.jl")
include("Query_test_yd.jl")
include("CS_generic.jl")
include("Test_degeneracy_yd.jl")

# Was 0 indexed, convert to 1-indexed.
# Remove any self-loops.

CS_DBLP_FOLDER = "../CaseStudy/DBLP/"
CS_DBLP_RAW_FOLDER = folderString(CS_DBLP_FOLDER, "Raw")
CS_DBLP_IN_FOLDER = folderString(CS_DBLP_FOLDER, "IN")

DBLP_NAME_FILE = "ent.author"
DBLP_CI_FILE = "csdblp.in"
DBLP_AUTHOR_TOTAL = 1824701
DBLP_COLLAB_TOTAL = 8344615

# TODO:
# Read raw to make an array of all nodes so that can index -> article title.

function LoadDBLPNameAsArray()
    io_read = open(string(CS_DBLP_RAW_FOLDER,DBLP_NAME_FILE))
    names = emptyStringArray(DBLP_AUTHOR_TOTAL)
    while !eof(io_read)
        line = readline(io_read)
        if !startswith(line, "%")
            lineSplit = split(line, " ")
            ind, name = parse(Int64, lineSplit[1]), lineSplit[2]
            names[ind] = name
        end
    end
    close(io_read)
    return names
end

println("Reading DBLP citation data...")
B = readIN(DBLP_CI_FILE, CS_DBLP_IN_FOLDER)
P = toTransitionGraph(B)
BW = readMulti("out.dblp_coauthor", DBLP_AUTHOR_TOTAL, DBLP_COLLAB_TOTAL, CS_DBLP_RAW_FOLDER)
allNames = LoadDBLPNameAsArray()

# V = 95485
# C = GenerateUserInputSet(B,V,2,4)

function GetRefinedSetDBLP(C::Vector{Int64})
    # As an example, say V = 95485, for getting info from the raw citation data given this index (note -1 for the raw file):
    # grep -n "#index95484" DBLP-citation-Jan8.txt
    # Say you get row number = 6809987, then:
    # Raw$ awk 'FNR>=6809987 && FNR<=6810020' DBLP-citation-Jan8.txt
    return GetRefinedSet(B, C, allTitles)
end

