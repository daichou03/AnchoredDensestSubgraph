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
include("Test_yd.jl")
include("Query_test_yd.jl")

# Was 0 indexed, convert to 1-indexed.
# Remove any self-loops.

function ConvertDBLPCitationToIN(FileName::AbstractString, OutputFileName::AbstractString, RawDirectory::String="../CaseStudy/Raw/", OutputDirectory::String="../CaseStudy/IN/")
    io_read = open(string(RawDirectory,FileName))
    N = parse(Int64, readline(io_read))
    M = 0
    node1 = -1
    v1 = Int64[]
    v2 = Int64[]
    while !eof(io_read)
        line = readline(io_read)
        if startswith(line, "#index")
            node1 = parse(Int64, line[7:length(line)]) + 1
        elseif startswith(line, "#%")
            node2 = parse(Int64, line[3:length(line)]) + 1
            if node1 != node2
                M += 1
                push!(v1, node1)
                push!(v2, node2)
            end
        end
    end
    close(io_read)
    # Write
    mkpath(OutputDirectory)
    io_write = open(string(OutputDirectory,OutputFileName), "w")
    write(io_write, string(N, " ", M, "\n"))
    for i = 1:M
        write(io_write, string(v1[i], " ", v2[i], "\n"))
    end
    close(io_write)
end

# TODO:
# Read raw to make an array of all nodes so that can index -> article title.

function LoadDBLPTitleAsArray(FileName::AbstractString, RawDirectory::String="../CaseStudy/Raw/")
    io_read = open(string(RawDirectory,FileName))
    titles = String[]
    while !eof(io_read)
        line = readline(io_read)
        if startswith(line, "#*")
            push!(titles, chop(line,head=2,tail=0))
        end
    end
    close(io_read)
    return titles
end

# TODO:
# In case we find a small CC, report. make a function that checks if this node is in a CC of at least k nodes. Maybe not necessary if we always cherry pick starting nodes.

DBLP_RAW_FILE = "DBLP-citation-Jan8.txt"
DBLP_CI_FILE = "dblpciv4.in"

# ConvertDBLPCitationToIN("DBLP-citation-Jan8.txt", DBLP_CI_FILE)

println("Reading DBLP citation data...")
B = readIN(DBLP_CI_FILE)
allTitles = LoadDBLPTitleAsArray(DBLP_RAW_FILE)

C = GenerateUserInputSet(B,V,2,4)

# V = 95485
function GetRefinedSet(C::Vector{Int64})
    # As an example, say V = 95485, for getting info from the raw citation data given this index (note -1 for the raw file):
    # grep -n "#index95484" DBLP-citation-Jan8.txt
    # Say you get row number = 6809987, then:
    # Raw$ awk 'FNR>=6809987 && FNR<=6810020' DBLP-citation-Jan8.txt
    
    R = GenerateReferenceSetFixedWalks(B,C)
    refined = LocalAnchoredDensestSubgraph(B,R).source_nodes
    println("Articles from refined set: ")
    println("------------")
    for i in refined
        println(string(allTitles[i]))
    end
    println("------------")
    return refined
end