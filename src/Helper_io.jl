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

DIR_EXAMPLE_SCC = "../Example_SCC/"

function readRaw(FileName::AbstractString, N::Int, M::Int, Chance::Float64=1.0, Directory::String=DIR_EXAMPLE_SCC)
    f = open(string(Directory,FileName))
    return doReadIN(f, N, M, Chance)
end

function readRaw(FileName::AbstractString, N::Int, M::Int, Directory::String=DIR_EXAMPLE_SCC)
    return readRaw(FileName, N, M, 1.0, Directory)
end

function readIN(FileName::AbstractString, Chance::Float64=1.0, Directory::String=DIR_EXAMPLE_SCC)
    f = open(string(Directory,FileName))
    header = split(readline(f))
    n = parse(Int,header[1])
    m = parse(Int,header[2])
    return doReadIN(f, n, m, Chance, false)
end

function readIN(FileName::AbstractString, Directory::String=DIR_EXAMPLE_SCC)
    return readIN(FileName, 1.0, Directory)
end

COMMENT_HEADERS = ['#','%']

# Any header info has been read, File now only has body content remaining
function doReadIN(File::IOStream, N::Int64, M::Int64, Chance::Float64, Weighted::Bool)
    ei = zeros(Int64, M*2)
    ej = zeros(Int64, M*2)
    if Weighted
        weights = zeros(Float64, M*2)
    end
    edgeIndex = map(x->Dict(), 1:N)
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
        weight = Weighted ? ((length(line) >= 3) ? parse(Float64, line[3]) : 1.0) : 0.0
        if haskey(edgeIndex[v1], v2)
            if Weighted
                weights[edgeIndex[v1][v2]*2-1] += weight
                weights[edgeIndex[v1][v2]*2] += weight
            end
        else
            count += 1
            ei[count*2-1] = v1
            ej[count*2-1] = v2
            ei[count*2] = v2
            ej[count*2] = v1
            if Weighted
                weights[count*2-1] = weight
                weights[count*2] = weight
            end
            edgeIndex[v1][v2] = count
        end
    end
    close(File)
    A = sparse(ei[1:count*2], ej[1:count*2], (Weighted ? weights[1:count*2] : ones(Float64, count*2)), N, N)
    return A
end


# Undirected and assumes input is undirected
# N, M can be lowerbound rather than accurate.
function readMulti(FileName::AbstractString, N::Int64, M::Int64, Directory::String=DIR_EXAMPLE_SCC)
    io = open(string(Directory,FileName))
    return doReadIN(io, N, M, 1.0, true)
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
