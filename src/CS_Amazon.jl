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

##################################
# For simple effectiveness tests #
##################################

# I/O R and V

function ExportSimpleRs(Vs::Vector{Int64}, Rs::Any, Name::String="SimpleRs")
    folder = string(CS_AMAZON_FOLDER, Name, "/")
    mkpath(folder)
    io_v = open(string(folder,"v.txt"), "w")
    io_r = open(string(folder,"r.txt"), "w")
    for i in 1:length(Vs)
        write(io_v, string(Vs[i], "\n"))
        write(io_r, string(join(Rs[i], ","), "\n"))
    end
    close(io_v)
    close(io_r)
end

function ImportSimpleRs(Name::String="SimpleRs")
    folder = string(CS_AMAZON_FOLDER, Name, "/")
    vs = Int64[]
    rs = Any[]

    io_v = open(string(folder,"v.txt"))
    while !eof(io_v)
        push!(vs, parse(Int64, readline(io_v)))
    end
    close(io_v)

    io_r = open(string(folder,"r.txt"))
    while !eof(io_r)
        push!(rs, map(x->parse(Int64, x), split(readline(io_r), ",")))
    end
        
    println(string(length(vs), " Rs imported."))
    return vs, rs
end

# I/O result set and time

function ExportSimpleResults(Ss::Any, Times::Any, TestName::String, DataName::String, AlgName::String)
    folder_s = string(CS_AMAZON_FOLDER, "Result/", TestName, "/", DataName, "/")
    folder_time = string(CS_AMAZON_FOLDER, "Time/", TestName, "/", DataName, "/")
    mkpath(folder_s)
    mkpath(folder_time)
    io_s = open(string(folder_s, AlgName, ".txt"), "w")
    io_time = open(string(folder_time, AlgName, ".txt"), "w")
    for j in 1:length(Ss)
        write(io_s, string(join(Ss[j], ","), "\n"))
        write(io_time, string(Times[j], "\n"))
    end
    close(io_s)
    close(io_time)
end

function ImportSimpleResults(TestName::String, DataName::String, AlgName::String)
    ss = Any[]
    io_s = open(string(string(CS_AMAZON_FOLDER, "Result/", TestName, "/", DataName, "/"),AlgName,".txt"))
    while !eof(io_s)
        push!(ss, map(x->parse(Int64, x), split(readline(io_s), ",")))
    end
    close(io_s)

    times = Any[]
    io_time = open(string(string(CS_AMAZON_FOLDER, "Time/", TestName, "/", DataName, "/"),AlgName,".txt"))
    while !eof(io_time)
        push!(times, parse(Float64, readline(io_time)))
    end
    close(io_time)

    println(string(length(ss), " results imported."))
    return ss, times
end
