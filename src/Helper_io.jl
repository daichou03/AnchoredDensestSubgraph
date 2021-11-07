using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra

# read "IN" format.
# Based on MatrixNetworks::readSMAT.
# IN format assumes the graph is undirectional, unweighted, 1-indexed, has header for n and m.

# Example for IN format:
# 5 7
# 1 2
# 1 3
# 1 4
# 2 3
# 2 4
# 3 4
# 3 5
# 4 5

# Note: if reading data with multi-edges, M should count duplicate edges multiple times (not once).

DIR_EXAMPLE_SCC = "../Example_SCC/"

function readRaw(FileName::AbstractString, N::Int, M::Int, Chance::Float64=1.0, Directory::String=DIR_EXAMPLE_SCC)
    f = open(string(Directory,FileName))
    return doReadIN(f, N, M, Chance, false, false)
end

function readRaw(FileName::AbstractString, N::Int, M::Int, Directory::String=DIR_EXAMPLE_SCC)
    return readRaw(FileName, N, M, 1.0, Directory)
end

function readIN(FileName::AbstractString, Chance::Float64=1.0, Directory::String=DIR_EXAMPLE_SCC)
    f = open(string(Directory,FileName))
    header = split(readline(f))
    n = parse(Int,header[1])
    m = parse(Int,header[2])
    return doReadIN(f, n, m, Chance, false, false)
end

function readIN(FileName::AbstractString, Directory::String=DIR_EXAMPLE_SCC)
    return readIN(FileName, 1.0, Directory)
end

COMMENT_HEADERS = ['#','%']

# Any header info has been read, File now only has body content remaining
# Multi: The input might have multi edges, combine them into weights (Sparse function does this by default). If false, force all weights = 1 afterwards.
# Weighted: Attempts to read weight too. Note if Weighted = true, Multi is treated as true.
function doReadIN(File::IOStream, N::Int64, M::Int64, Chance::Float64, Multi::Bool, Weighted::Bool)
    ei = zeros(Int64, M*2)
    ej = zeros(Int64, M*2)
    weights = []
    if Weighted
        weights = zeros(Float64, M*2)
    end
    count = 0
    while !eof(File)
        lineRaw = readline(File)
        if length(lineRaw) == 0 || (lineRaw[1] in COMMENT_HEADERS) || (Chance < rand())
            continue
        end
        lineRaw = replace(lineRaw, "\t"=>" ")
        line = split(lineRaw, " ")
        v1 = parse(Int, line[1])
        v2 = parse(Int, line[2])
        if v1 == v2
            continue
        elseif v1 > v2
            v1, v2 = v2, v1
        end
        weight = (Weighted && length(line) >= 3) ? parse(Float64, line[3]) : 1.0

        count += 1
        ei[count*2-1] = v1
        ej[count*2-1] = v2
        ei[count*2] = v2
        ej[count*2] = v1
        if Weighted
            weights[count*2-1] = weight
            weights[count*2] = weight
        end
    end
    close(File)
    A = sparse(ei[1:count*2], ej[1:count*2], (Weighted ? weights[1:count*2] : ones(Float64, count*2)), N, N)
    if !Weighted && !Multi
        # Force weight = 1
        for i = 1:length(A.nzval)
            A.nzval[i] = 1
        end
    end
    return A
end


# Assumes input is undirected and multi edge, try to merge these edges into a weighted graph.
# N, M can be lowerbound rather than accurate.
function readMulti(FileName::AbstractString, N::Int64, M::Int64, Directory::String=DIR_EXAMPLE_SCC)
    io = open(string(Directory,FileName))
    return doReadIN(io, N, M, 1.0, true, false)
end

function exportIN(B::SparseMatrixCSC, FileName::String, Directory::String="../Example_SCC/")
    io = open(string(Directory,FileName), "w")
    N = size(B,1)
    write(io, string(size(B,1)," ",Int64(nnz(B)/2),"\n"))
    for i = 1:N
        indices = B[:,i].nzind
        indices = indices[searchsortedfirst(indices, i) : length(indices)]
        for j = 1:length(indices)
            write(io, string(i," ",indices[j],"\n"))
        end
    end
    close(io)
end

# Generate multiple graphs, each one has half edge as the previous.
function ExportHalfEdgeGraphs(GraphName::String, Iteration::Integer=5)
    println(string("Graph: ",GraphName))
    iter = 0
    while iter < Iteration
        iter += 1
        println(string("Iteration: ", iter))
        g = RetrieveLargestConnectedComponent(readIN(string(GraphName, ".in"), 0.5^iter))
        # println(string("Subgraph size = ", size(g, 1), ", now exporting..."))
        exportIN(g, string(GraphName, "-H", iter, ".in"))
    end
end

function BulkExportHalfEdgeGraphs(dataset_names::Array{String,1}, Iteration::Integer=5)
    for ds_name in dataset_names
        ExportHalfEdgeGraphs(ds_name, Iteration)
    end
end

#"lastfm","deezer","orkut","livejournal","dblp","youtube","amazon","github","astroph","condmat","grqc","hepph","hepth","brightkite","catster","hamster","douban","gowalla","douban","gowalla","gowalla","douban","gowalla"


# In general, for new data, load it, take its largest connected component, and then output it back to /Example_SCC/ for later use.
# Example:

# epinion = RetrieveLargestConnectedComponent(readIN("soc-Epinions1.in"), "../Example/")
# exportIN(epinion, "epinion.in")

# exportIN(RetrieveLargestConnectedComponent(readIN("gowalla_edges.in", "../Example_preprocessed/")), "gowalla.in")
