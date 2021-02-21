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

function readIN(FileName::AbstractString, Directory::String="../Example_SCC/")
    f = open(string(Directory,FileName))
    header = readline(f)
    headerparts = split(header)
    nedges = parse(Int,headerparts[2])
    ei = zeros(Int64,nedges*2)
    ej = zeros(Int64, nedges*2)
    @inbounds for i = 1:nedges
        curline = readline(f)
        parts = split(curline)
        ei[2*i-1] = parse(Int, parts[1])
        ej[2*i-1] = parse(Int, parts[2])
        ei[2*i] = parse(Int, parts[2])
        ej[2*i] = parse(Int, parts[1])
    end
    close(f)
    A = sparse(ei, ej, ones(Float64, nedges*2),
               parse(Int,headerparts[1]), 
               parse(Int,headerparts[1]))
    return A
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

# In general, for new data, load it, take its largest connected component, and then output it back to /Example_SCC/ for later use.
# Example:

# epinion = RetrieveLargestConnectedComponent(readIN("soc-Epinions1.in"), "../Example/")
# exportIN(epinion, "epinion.in")

# exportIN(RetrieveLargestConnectedComponent(readIN("com-orkut.ungraph.in", "../Example/"), "../Example/"), "orkut.in")
