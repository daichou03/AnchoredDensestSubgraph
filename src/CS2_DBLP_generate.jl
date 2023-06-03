using CSV
using DataFrames
using StatsBase
include("Utils.jl")
include("Utils_io.jl")
include("LP_consts.jl")
include("Utils_graph.jl")
include("CS_DBLP.jl")
include("CS2_generic.jl")
include("LP_algorithm.jl")

dataName = "csdblp"
FOLDER_CS_DBLP_CANDIDATE_LP = folderString(FOLDER_CS_DBLP, "candidate", "LP")

V_JW = 16028
N_JW = GetAdjacency(B, V_JW, true)
N2R_JW = GetComponentAdjacency(B, N_JW, false)

# Very specific - make use of results produced by CS_DBLP_LA.CandidateSearchStore to produce R.
function collectCandidateAsR()
    Vs, Rs = [], []
    files = readdir(FOLDER_CS_DBLP_CANDIDATE_LP)
    
    for file in files
        if all(isdigit, split(file, ".")[1])
            file_path = joinpath(FOLDER_CS_DBLP_CANDIDATE_LP, file)
            io = open(file_path, "r")
            push!(Vs, parse(Int64, readline(io)))
            push!(Rs, parse.(Int64, split(readline(io), ",")))
            close(io)
        end
    end
    writeAnchors(dataName, "baseline", Rs)
    # Write Vs to the same folder
    file = open(string(FOLDER_CS_DBLP_CANDIDATE_LP, "Vs.txt"), "w")
    for V in Vs
        write(io, string(V, "\n"))
    end
    close(file)
end
