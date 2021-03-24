using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra

# ydai992: read "IN" format.
# Based on MatrixNetworks::readSMAT.
# IN format assumes the graph is undirectional, unweighted, 1-indexed, has header for m and n.

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

function readIN(FileName::AbstractString, Chance::Float64=1.0, Directory::String="../Example_SCC/")
    f = open(string(Directory,FileName))
    header = split(readline(f))
    nedges = parse(Int,header[2])
    ei = zeros(Int64, nedges*2)
    ej = zeros(Int64, nedges*2)
    count = 0
    @inbounds for i = 1:nedges
        curline = readline(f)
        if Chance >= 1 || Chance >= rand()
            count += 1
            parts = split(curline)
            ei[2*count-1] = parse(Int, parts[1])
            ej[2*count-1] = parse(Int, parts[2])
            ei[2*count] = parse(Int, parts[2])
            ej[2*count] = parse(Int, parts[1])
        end
    end
    close(f)
    A = sparse(ei[1:count*2], ej[1:count*2], ones(Float64, count*2),
               parse(Int,header[1]), 
               parse(Int,header[1]))
    return A
end

function readIN(FileName::AbstractString, Directory::String="../Example_SCC/")
    return readIN(FileName, 1.0, Directory)
end

function exportIN(B::SparseMatrixCSC, FileName::String, Directory::String="../Example_SCC/")
    io = open(string(Directory,FileName), "w")
    N = size(B,1)
    write(io, string(size(B,1)," ",Int64(nnz(B)/2),"\n"))
    for i = 1:N
        indices = B[i,:].nzind
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

# Was 0 indexed, convert to 1-indexed.
# Remove any self-loops.

function ConvertDBLPCitationToIN(FileName::AbstractString, RawDirectory::String="../CaseStudy/Raw/", OutputDirectory::String="../CaseStudy/IN/")
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
        else if startswith(line, "#%")
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
    io_write = open(string(OutputDirectory,FileName), "w")
    write(io_write, string(N, " ", M, "\n"))
    for i = 1:M
        write(io_write, v1[i], " ", v2[i], "\n")
    end
    close(io_write)
end

#"lastfm","deezer","orkut","livejournal","dblp","youtube","amazon","github","astroph","condmat","grqc","hepph","hepth","brightkite","catster","hamster","douban","gowalla","douban","gowalla","gowalla","douban","gowalla"


# In general, for new data, load it, take its largest connected component, and then output it back to /Example_SCC/ for later use.
# Example:

# epinion = RetrieveLargestConnectedComponent(readIN("soc-Epinions1.in"), "../Example/")
# exportIN(epinion, "epinion.in")

# exportIN(RetrieveLargestConnectedComponent(readIN("gowalla_edges.in", "../Example_preprocessed/")), "gowalla.in")
