using SparseArrays
using MAT
using MatrixNetworks
using LinearAlgebra
using StatsBase
using Random
using Base
include("Helper_io.jl")
include("Utils.jl")
include("Graph_utils_yd.jl")

# This part of the code should be made safe to be included by both our algorithm and any competitor algorithms.

CS_AMAZON_FOLDER = "../CaseStudy/Amazon/"

AMAZON_GRAPH_FILE = string(CS_AMAZON_FOLDER, "IN/com-amazon.ungraph.in")

println("Reading Amazon data...")
B = readIN(AMAZON_GRAPH_FILE)
P = toTransitionGraph(B)

function ExportIndicesByDegree(Last::Int64=40)
    ios = Any[]
    for i = 1:(Last+1)
        append!(ios, 0)
        ios[i] = open(string(CS_AMAZON_FOLDER,string("Helper-ind/",i,".txt")), "w")
    end
    for v in 1:size(B,1)
        deg = GetDegree(B, v)
        if deg > Last
            write(ios[Last+1], string(v, "\n"))
        elseif deg > 0
            write(ios[deg], string(v, "\n"))
        end
    end
    for i = 1:(Last+1)
        close(ios[i])
    end
end

function ReadIndicesByDegree(Last::Int64=40)
    inds = Any[]
    for i = 1:(Last+1)
        io = open(string(CS_AMAZON_FOLDER,string("Helper-ind/",i,".txt")))
        append!(inds, 0)
        inds[i] = []
        while !eof(io)
            append!(inds[i], parse(Int64, readline(io)))
        end
        close(io)
    end
    return inds
end

# For Stratified effectiveness tests

function ExportRs(Rs::Any, Name::String="0")
    folder = string(CS_AMAZON_FOLDER, "R-", Name, "/")
    mkpath(folder)
    for i in 1:length(Rs)
        io = open(string(folder,i,".txt"), "w")
        for j in 1:length(Rs[i])
            write(io, string(join(Rs[i][j], ","), "\n"))
        end
        close(io)
    end
end

function ImportRs(Name::String="0")
    folder = string(CS_AMAZON_FOLDER, "R-", Name, "/")
    rs = Any[]
    i = 1
    filename = string(folder,i,".txt")
    while isfile(filename)
        io = open(filename)
        append!(rs, 0)
        rs[i] = []
        j = 1
        while !eof(io)
            append!(rs[i], 0)
            rs[i][j] = map(x->parse(Int64, x), split(readline(io), ","))
            j += 1
        end
        close(io)
        i += 1
        filename = string(folder,i,".txt")
    end
    println(string(i, " set of Rs imported."))
    return rs
end
